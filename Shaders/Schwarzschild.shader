Shader "Unlit/Schwarzschild"
{
    Properties
    {
        [HideInInspector] [Toggle(ATTACHED_TO_CUBE)] _AttachedToCube ("Attached To Cube", Float) = 0
        [Header(Configuration)]
        _QuadScale ("Quad Scale", Float) = 3.0
        [HideInInspector] _CubeScale ("Cube Scale", Float) = 100.0

        [Header(Spacetime)]
        _c ("Speed Of Light", Float) = 1.0

        [Header(Blackhole)]
        _a ("Schwarzschild Radius", Float) = 1.0

        [Header(Calculation Parameters)]
        _N ("N Steps", Int) = 500
        _dtau ("Step Size", Float) = 0.333
        _EscapeVelocity ("Escape Velocity", Range(0, 1)) = 0.8
        _MaxWinding ("Max Winding Of Light", Float) = 3.0

        [Header(Ring)]
        _RingRadius ("Ring Radius", Float) = 8.0
        [NoScaleOffset] _RedShiftTex ("Red Shift Texture", 2D) = "white" {}
        [Header(Ring Noise)]
        _RingNoise_r ("N Divisions Of r", Int) = 8
        _RingNoise_phi ("N Divisions Of phi", Int) = 3

        [Header(Bloom)]
        _sigma ("sigma", Float) = 1.25
        _StepWidth ("Step Width", Float) = 8.0
        _Threshold ("Threshold", Float) = 0.8
        _Suppression ("Suppression", Float) = 0.7
    }

    SubShader
    {
        LOD 100

        CGINCLUDE
            #pragma vertex vert

            #pragma shader_feature DRAW_STARS

            #include "UnityCG.cginc"
            #include "Functions.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
            };

            struct vertout
            {
                float3 ray_pos : POSITION1;
                float4 vertex : SV_POSITION;
                float4 uv : TEXCOORD0;
            };

            float _QuadScale, _CubeScale;
            float _RingRadius;

            float _sigma;
            float _StepWidth;
            float _Threshold;
            float _Suppression;

            vertout vert (appdata v)
            {
                vertout o;

                #if ATTACHED_TO_CUBE
                    o.vertex = float4(- _CubeScale * v.vertex.xyz, 1);
                    o.ray_pos = zinv(o.vertex.xyz);
                    o.vertex = UnityObjectToClipPos(o.vertex);
                #else
                    float3 obj = UNITY_MATRIX_T_MV[3].xyz;
                    o.vertex = float4(obj + _QuadScale * _RingRadius * v.vertex, 1);
                    o.ray_pos = zinv(mul(o.vertex, UNITY_MATRIX_IT_MV).xyz);
                    o.vertex = mul(UNITY_MATRIX_P, o.vertex);
                #endif

                o.uv = ComputeGrabScreenPos(o.vertex);
                return o;
            }
        ENDCG

        Tags { "RenderType"="Transparent" }
        Blend SrcAlpha OneMinusSrcAlpha

        Pass
        {
            CGPROGRAM
            #pragma fragment frag

            struct fragin
            {
                float3 ray_pos : POSITION1;
                float4 uv : TEXCOORD0;
            };

            struct fragout
            {
                float4 color : SV_Target;
                float depth : SV_Depth;
            };

            float _c;
            float _a;
            int _N;
            float _dtau;
            float _EscapeVelocity;
            float _MaxWinding;

            sampler2D _RedShiftTex;
            int _RingNoise_r, _RingNoise_phi;

            sampler2D _CameraDepthTexture;

            float g00 (float r)
            {
                return 1 - _a / r;
            }

            float delw (float c0, float r)
            {
                return - c0 / g00(r);
            }

            float del2r (float c0, float omega0, float r)
            {
                return (r - 1.5 * _a) * pow(omega0 / (r * r), 2);
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

                float u0w = _c / rtg00r0;
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
                return _dtau * r_ur(rur[ur_], del2r(c0, omega0, rur[r_]));
            }

            r_ur runge_kutta (float c0, float omega0, r_ur rur)
            {
                float half_dtau = 0.5 * _dtau;

                r_ur k1 = half_dtau * r_ur(rur[ur_],           del2r(c0, omega0, rur[r_]));
                r_ur k2 = half_dtau * r_ur(rur[ur_] + k1[ur_], del2r(c0, omega0, rur[r_] + k1[r_]));
                r_ur k3 =     _dtau * r_ur(rur[ur_] + k2[ur_], del2r(c0, omega0, rur[r_] + k2[r_]));
                r_ur k4 = half_dtau * r_ur(rur[ur_] + k3[ur_], del2r(c0, omega0, rur[r_] + k3[r_]));

                return (1.0 / 3.0) * (k1 + 2*k2 + k3 + k4);
            }

            void init_integral (float c0, float omega0, r_ur rur, out w_phi prev)
            {
                prev = _dtau * w_phi(delw(c0, rur[r_]), delphi(omega0, rur[r_]));
            }

            w_phi trapezoid (float c0, float omega0, r_ur rur, w_phi wphi, inout w_phi prev)
            {
                w_phi now = _dtau * w_phi(delw(c0, rur[r_]), delphi(omega0, rur[r_]));
                w_phi ret = 0.5 * (prev + now);
                prev = now;
                return ret;
            }

            w_phi simpson (int n, float c0, float omega0, r_ur rur, w_phi wphi, inout w_phi prev1, inout w_phi prev2)
            {
                w_phi now = _dtau * w_phi(delw(c0, rur[r_]), delphi(omega0, rur[r_]));
                w_phi ret = fmod(n, 2) ? 0.5 * (prev1 + now) : (1.0 / 3.0) * (- 0.5 * prev2 + 2.5 * prev1 + now);
                prev2 = prev1;
                prev1 = now;
                return ret;
            }

            fragout frag (fragin i)
            {
                fragout o;
                o.color = 0;
                o.depth = 1;

                float3 obs3 = zinv(mul(unity_WorldToObject, _WorldSpaceCameraPos));
                float3 u03_obs = _c * normalize(i.ray_pos - obs3);

                float3x3 rot = calc_rotation(obs3, u03_obs);
                float2 obs = mul(rot, obs3).xy;
                float2 u0_obs = mul(rot, u03_obs).xy;
                float2 proj_ey = mul(rot, float3(0, 1, 0)).xy;
                float eyphi = fmod(atan2(proj_ey.y, proj_ey.x) + 1.5 * PI, PI);
                float phi_shift = eyphi < 0.5 * PI ? 0.5 * PI : 0;
                float4x4 clip_inv_z_rot = mul(UNITY_MATRIX_VP, mul(unity_ObjectToWorld, f3x3to4x4(mul(ZINV_MATRIX, transpose(rot)))));

                float c0, omega0;
                r_ur rur;
                w_phi wphi;
                float rtg00r0;
                init(obs, u0_obs, c0, omega0, wphi[w_], rur[r_], rur[ur_], wphi[phi_], rtg00r0);
                float r0 = rur[r_];

                bool decided = false;
                w_phi dwphi_prev1, dwphi_prev2;
                r_ur drur;
                w_phi dwphi;
                float diff_phi, diff_phi_prev;
                float r_lerp, r_prev;
                w_phi wphi_lerp, wphi_prev;
                float s;
                float3 ring_pos;
                bool flag;
                fixed4 color;
                float max_depth = 0;
                float prev_depth = 0;
                float buffer_depth = UNITY_SAMPLE_DEPTH(tex2D(_CameraDepthTexture, f4to3(i.uv).xy));

                init_integral(c0, omega0, rur, dwphi_prev1);
                for (int n = 1; n <= _N & ! decided; n++) {
                    max_depth = max(prev_depth, max_depth);
                    r_prev = rur[r_];
                    wphi_prev = wphi;

                    //drur = euler(c0, omega0, rur);
                    drur = runge_kutta(c0, omega0, rur);
                    rur += drur;
                    //dwphi = trapezoid(c0, omega0, rur, wphi, dwphi_prev1);
                    dwphi = simpson(n, c0, omega0, rur, wphi, dwphi_prev1, dwphi_prev2);
                    wphi += dwphi;

                    flag = wphi[phi_] > 2 * PI * _MaxWinding;
                    o.color.rgb = flag ? o.color.a * o.color.rgb : o.color.rgb;
                    o.color.a = flag ? 1 : o.color.a;
                    o.depth = flag ? max_depth : o.depth;
                    decided = decided || flag;

                    decided = decided || rur[ur_] / rtg00r0 > _EscapeVelocity * _c;

                    s = (_a - r_prev) / drur[r_];
                    wphi_lerp[phi_] = wphi_prev[phi_] + dwphi[phi_] * s;

                    flag = rur[r_] < _a;
                    o.color.rgb = flag ? o.color.a * o.color.rgb : o.color.rgb;
                    o.color.a = flag ? 1 : o.color.a;
                    o.depth = flag ? max(clip2depth(mul(clip_inv_z_rot, float4(r_lerp * float2(cos(wphi_lerp[phi_]), sin(wphi_lerp[phi_])), 0, 1))), max_depth) : o.depth;
                    decided = decided || flag;

                    diff_phi_prev = fmod(wphi_prev[phi_] + phi_shift, PI) - phi_shift - eyphi;
                    diff_phi = diff_phi_prev + dwphi[phi_];
                    s = - diff_phi_prev / dwphi[phi_];
                    r_lerp = r_prev + drur[r_] * s;
                    wphi_lerp = wphi_prev + dwphi * s;
                    ring_pos = r_lerp * zinv(mul(float3(cos(wphi_lerp[phi_]), sin(wphi_lerp[phi_]), 0), rot));

                    flag =
                        diff_phi * diff_phi_prev < 0 &&
                        3 * _a < r_lerp && r_lerp < _RingRadius;

                    color.rgb = tex2Dlod(_RedShiftTex, float4(3 * _a / r_lerp, _a / r0, 0, 0)).rgb;
                    color.a = saturate(perlin(float2(_RingRadius - 3 * _a, 2 * PI), int2(_RingNoise_r, _RingNoise_phi), float2(r_lerp - 3 * _a, atan2(ring_pos.x, ring_pos.z) + PI + _Time.z)) * 2 + 0.8);
                    color.a *= smoothstep(3 * _a, 3 * _a + 0.3, r_lerp) * smoothstep(_RingRadius, _RingRadius - 0.3, r_lerp);

                    o.color.rgb = flag ? o.color.a * o.color.rgb + (1 - o.color.a) * color.rgb : o.color.rgb;
                    o.color.a = flag ? 1 - (1 - o.color.a) * (1 - color.a) : o.color.a;
                    o.depth = flag ? max(clip2depth(UnityObjectToClipPos(ring_pos)), max_depth) : o.depth;
                    decided = decided || (flag && o.color.a > 0.99);

                    prev_depth = clip2depth(mul(clip_inv_z_rot, float4(rur[r_] * float2(cos(wphi[phi_]), sin(wphi[phi_])), 0, 1)));
                    decided = decided || prev_depth > buffer_depth;
                }
                flag = rur[r_] < 1.5 * _a;
                o.color.rgb = flag ? o.color.a * o.color.rgb : o.color.rgb;
                o.color.a = flag ? 1 : o.color.a;
                o.depth = flag ? max_depth : o.depth;

                clip(o.color.a - 0.01);
                return o;
            }
            ENDCG
        }

        Tags { "RenderType"="Transparent" }
        Blend SrcAlpha OneMinusSrcAlpha
        ZTest Always
        ZWrite Off

        GrabPass {}

        Pass
        {
            CGPROGRAM
            #pragma fragment frag

            struct fragin
            {
                float4 uv : TEXCOORD0;
            };

            sampler2D _GrabTexture;
            float4 _GrabTexture_TexelSize;

            fixed4 frag (fragin i) : SV_Target
            {
                float2 uv = f4to3(i.uv).xy;
                fixed4 color = tex2D(_GrabTexture, uv);
                color.xyz += xblur(_GrabTexture, gaussian_filter5(_sigma), uv, _StepWidth * _GrabTexture_TexelSize.x, _Threshold, _Suppression);
                color.a = 1;
                return color;
            }
            ENDCG
        }

        GrabPass {}

        Pass
        {
            CGPROGRAM
            #pragma fragment frag

            struct fragin
            {
                float4 uv : TEXCOORD0;
            };

            sampler2D _GrabTexture;
            float4 _GrabTexture_TexelSize;

            fixed4 frag (fragin i) : SV_Target
            {
                float2 uv = f4to3(i.uv).xy;
                fixed4 color = tex2D(_GrabTexture, uv);
                color.xyz += yblur(_GrabTexture, gaussian_filter5(_sigma), uv, _StepWidth * _GrabTexture_TexelSize.y, _Threshold, _Suppression);
                color.a = 1;
                return color;
            }
            ENDCG
        }
    }
}
