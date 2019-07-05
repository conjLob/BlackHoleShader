Shader "Unlit/Schwarzschild"
{
    Properties
    {
        [NoScaleOffset]
        _RedShiftTex ("Red Shift Texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags { "RenderType"="Transparent" }
        Blend SrcAlpha OneMinusSrcAlpha
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "UnityCG.cginc"
            #include "Functions.cginc"

            #define a 1.0
            #define c 1.0
            #define dtau 0.333
            #define N 500
            #define ESCAPE_VELOCITY 0.8
            #define MAX_WINDING 3
            #define RING_SCALE 8.0

            struct appdata
            {
                float4 vertex : POSITION;
            };

            struct v2f
            {
                float3 raypos : POSITION1;
                float4 vertex : SV_POSITION;
            };

            sampler2D _RedShiftTex;

            v2f vert (appdata v)
            {
                v2f o;
                float3 obj = UNITY_MATRIX_T_MV[3].xyz;
                o.vertex = float4(obj + RING_SCALE / 5.0 * 1.5 * float3(v.vertex.xz, 0), 1);

                o.raypos = zinv(mul(o.vertex, UNITY_MATRIX_IT_MV).xyz);
                o.vertex = mul(UNITY_MATRIX_P, o.vertex);
                return o;
            }

            float g00 (float r)
            {
                return 1 - a / r;
            }

            float delw (float c0, float r)
            {
                return - c0 / g00(r);
            }

            float del2r (float c0, float omega0, float r)
            {
                return (r - 1.5 * a) * pow(omega0 / (r * r), 2);
            }

            float delphi (float omega0, float r)
            {
                return omega0 / (r * r);
            }

            float3x3 calc_rotation (float3 obs, float3 u0_obs)
            {
                float3 er = float3(normalize(obs.xy - obs.z / u0_obs.z * u0_obs.xy), 0);
                float3 ephi = float3(-er.y, er.x, 0);
                float3 proj_u0_obs = normalize(float3(dot(u0_obs, ephi) * ephi.xy, u0_obs.z));
                float cos_ = dot(proj_u0_obs, ephi);
                float sin_ = dot(cross(ephi, proj_u0_obs), er);

                return rotation_tri(er, cos_, -sin_);
            }

            void init (
                float2 obs, float2 u0_obs,
                out float c0, out float omega0,
                out float w0, out float r0, out float u0r, out float phi0,
                out float rtg00r0
            ) {
                w0 = 0;
                r0 = length(obs);
                phi0 = fmod(atan2(obs.y, obs.x) + 2 * PI, 2 * PI);

                float g00r0 = g00(r0);
                rtg00r0 = sqrt(g00r0);
                float2 er_obs = obs / r0;
                float2 ephi_obs = float2(-er_obs.y, er_obs.x);

                float u0w = c / rtg00r0;
                u0r = rtg00r0 * dot(u0_obs, er_obs);
                float u0phi = dot(u0_obs, ephi_obs) / r0;

                c0 = - g00r0 * u0w;
                omega0 = r0 * r0 * u0phi;
            }

            typedef float2 r_ur;
            #define r_  0
            #define ur_ 1
            typedef float2 w_phi;
            #define w_   0
            #define phi_ 1

            r_ur euler (float c0, float omega0, r_ur rur)
            {
                return dtau * r_ur(rur[ur_], del2r(c0, omega0, rur[r_]));
            }

            r_ur runge_kutta (float c0, float omega0, r_ur rur)
            {
                float half_dtau = 0.5 * dtau;

                r_ur k1 = half_dtau * r_ur(rur[ur_],           del2r(c0, omega0, rur[r_]));
                r_ur k2 = half_dtau * r_ur(rur[ur_] + k1[ur_], del2r(c0, omega0, rur[r_] + k1[r_]));
                r_ur k3 =      dtau * r_ur(rur[ur_] + k2[ur_], del2r(c0, omega0, rur[r_] + k2[r_]));
                r_ur k4 = half_dtau * r_ur(rur[ur_] + k3[ur_], del2r(c0, omega0, rur[r_] + k3[r_]));

                return (1.0 / 3.0) * (k1 + 2*k2 + k3 + k4);
            }

            void init_integral (float c0, float omega0, r_ur rur, out w_phi prev)
            {
                prev = dtau * w_phi(delw(c0, rur[r_]), delphi(omega0, rur[r_]));
            }

            w_phi trapezoid (float c0, float omega0, r_ur rur, w_phi wphi, inout w_phi prev)
            {
                w_phi now = dtau * w_phi(delw(c0, rur[r_]), delphi(omega0, rur[r_]));
                w_phi ret = 0.5 * (prev + now);
                prev = now;
                return ret;
            }

            w_phi simpson (int n, float c0, float omega0, r_ur rur, w_phi wphi, inout w_phi prev1, inout w_phi prev2)
            {
                w_phi now = dtau * w_phi(delw(c0, rur[r_]), delphi(omega0, rur[r_]));
                w_phi ret = fmod(n, 2) ? 0.5 * (prev1 + now) : (1.0 / 3.0) * (- 0.5 * prev2 + 2.5 * prev1 + now);
                prev2 = prev1;
                prev1 = now;
                return ret;
            }

            fixed4 frag (v2f input) : SV_Target
            {
                fixed4 color = 0;

                float3 obs3 = zinv(mul(unity_WorldToObject, _WorldSpaceCameraPos));
                float3 u03_obs = c * normalize(input.raypos - obs3);

                float3x3 rot = calc_rotation(obs3, u03_obs);
                float2 obs = mul(rot, obs3).xy;
                float2 u0_obs = mul(rot, u03_obs).xy;
                float2 proj_ey = mul(rot, float3(0, 1, 0)).xy;
                float eyphi = fmod(atan2(proj_ey.y, proj_ey.x) + 1.5 * PI, PI);

                float c0, omega0;
                r_ur rur;
                w_phi wphi;
                float rtg00r0;
                init(obs, u0_obs, c0, omega0, wphi[w_], rur[r_], rur[ur_], wphi[phi_], rtg00r0);
                float r0 = rur[r_];

                bool frag;
                bool decided = false;

                r_ur drur;
                w_phi dwphi;
                float phi_shift, diff_phi, diff_phi_prev;
                float r_lerp, r_prev;
                w_phi wphi_lerp, wphi_prev;
                float s;
                float2 uv;

                w_phi dwphi_prev1, dwphi_prev2;
                init_integral(c0, omega0, rur, dwphi_prev1);
                for (int n = 1; n <= N & ! decided; n++) {
                    r_prev = rur[r_];
                    wphi_prev = wphi;

                    //drur = euler(c0, omega0, rur);
                    drur = runge_kutta(c0, omega0, rur);
                    rur += drur;
                    //dwphi = trapezoid(c0, omega0, rur, wphi, dwphi_prev1);
                    dwphi = simpson(n, c0, omega0, rur, wphi, dwphi_prev1, dwphi_prev2);
                    wphi += dwphi;

                    phi_shift = 0.5 * PI * step(eyphi, 0.5 * PI);
                    diff_phi_prev = fmod(wphi_prev[phi_] + phi_shift, PI) - phi_shift - eyphi;
                    diff_phi = diff_phi_prev + dwphi[phi_];
                    s = - diff_phi_prev / dwphi[phi_];
                    r_lerp = r_prev + drur[r_] * s;
                    wphi_lerp = wphi_prev + dwphi * s;

                    //uv = (1.0 / RING_SCALE) * r_lerp * zinv(mul(float3(cos(wphi_lerp[phi_]), sin(wphi_lerp[phi_]), 0), rot)).xz + 0.5;
                    uv = mul(float3(cos(wphi_lerp[phi_]), sin(wphi_lerp[phi_]), 0), rot).xz;
                    color = tex2Dlod(_RedShiftTex, float4(3 * a / r_lerp, a / r0, 0, 0));
                    color.w = perlin(float2(3, 2 * PI), float2(5, 3), float2(r_lerp - 3 * a, atan2(uv.y, uv.x) + PI + _Time.z)) * 2 + 0.8;
                    color.w = color.w > 0.5;

                    decided = diff_phi * diff_phi_prev < 0 && 3 * a < r_lerp && r_lerp < RING_SCALE && color.w;
                    color *= decided;

                    frag = rur[r_] < a;
                    color.w = frag ? 1 : color.w;
                    decided = decided || frag;

                    frag = wphi[phi_] > 2 * PI * MAX_WINDING;
                    color.w = frag ? 1 : color.w;
                    decided = decided || frag;

                    decided = decided || rur[ur_] / rtg00r0 > ESCAPE_VELOCITY * c;
                }
                color.w = rur[r_] < 1.5 * a ? 1 : color.w;
                clip(color.w - 0.1);
                return color;
            }
            ENDCG
        }
    }
}
