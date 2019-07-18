#ifndef FUNCTIONS_INCLUDED
#define FUNCTIONS_INCLUDED

#define PI 3.14159265359

#define RANDOM_SEED_0 38.4382
#define RANDOM_SEED_1 43758.5453
#define RANDOM_SEED_2 8.4678
#define FLOAT4_SEED float4(12.9898, 78.233, 59.381, 36.697)

float random (float seed)
{
    return frac(RANDOM_SEED_1 * sin(seed + RANDOM_SEED_0));
}

#define RANDOM(seed, axis) random(dot(FLOAT4_SEED.axis, (seed).axis))

float2 random2 (float2 array)
{
    return float2(
        RANDOM(array, xy),
        RANDOM(array + float2(0, RANDOM_SEED_2), xy)
    );
}

float3 random3 (float3 array)
{
    return float3(
        RANDOM(array, xyz),
        RANDOM(array + float3(0, RANDOM_SEED_2, 0), xyz),
        RANDOM(array + float3(0, 0, RANDOM_SEED_2), xyz)
    );
}

float perlin (float2 size, int2 n, float2 uv)
{
    uv = frac(uv / size) * n;
    float2 p00 = floor(uv);
    float2 p11 = fmod(p00 + 1, n);
    float2 uv_ = smoothstep(0, 1, frac(uv));
    float2 g00 = random2(p00) * 2 - 1;
    float2 g10 = random2(float2(p11.x, p00.y)) * 2 - 1;
    float2 g01 = random2(float2(p00.x, p11.y)) * 2 - 1;
    float2 g11 = random2(p11) * 2 - 1;
    float2 v00 = uv_;
    float2 v10 = uv_ - float2(1, 0);
    float2 v01 = uv_ - float2(0, 1);
    float2 v11 = uv_ - float2(1, 1);
    return lerp(
        lerp(dot(g00, v00), dot(g10, v10), uv_.x),
        lerp(dot(g01, v01), dot(g11, v11), uv_.x),
        uv_.y
    );
}

float perlin (float3 size, int3 n, float3 pos)
{
    pos = frac(pos / size) * n;
    float3 p000 = floor(pos);
    float3 p111 = fmod(p000 + 1, n);
    float3 pos_ = smoothstep(0, 1, frac(pos));
    float3 g000 = random3(p000) * 2 - 1;
    float3 g100 = random3(float3(p111.x, p000.y, p000.z)) * 2 - 1;
    float3 g010 = random3(float3(p000.x, p111.y, p000.z)) * 2 - 1;
    float3 g001 = random3(float3(p000.x, p000.y, p111.z)) * 2 - 1;
    float3 g011 = random3(float3(p000.x, p111.y, p111.z)) * 2 - 1;
    float3 g101 = random3(float3(p111.x, p000.y, p111.z)) * 2 - 1;
    float3 g110 = random3(float3(p111.x, p111.y, p000.z)) * 2 - 1;
    float3 g111 = random3(p111) * 2 - 1;
    float3 v000 = pos_;
    float3 v100 = pos_ - float3(1, 0, 0);
    float3 v010 = pos_ - float3(0, 1, 0);
    float3 v001 = pos_ - float3(0, 0, 1);
    float3 v011 = pos_ - float3(0, 1, 1);
    float3 v101 = pos_ - float3(1, 0, 1);
    float3 v110 = pos_ - float3(0, 1, 1);
    float3 v111 = pos_ - float3(1, 1, 1);
    return lerp(
        lerp(
            lerp(dot(g000, v000), dot(g100, v100), pos_.x),
            lerp(dot(g010, v010), dot(g110, v110), pos_.x),
            pos_.y
        ),
        lerp(
            lerp(dot(g001, v001), dot(g101, v101), pos_.x),
            lerp(dot(g011, v011), dot(g111, v111), pos_.x),
            pos_.y
        ),
        pos_.z
    );
}

float square (float l, float r, float x)
{
    return step(l, x) * step(x, r);
}

float3 f4to3 (float4 v)
{
    return v.xyz / v.w;
}

float4x4 f3x3to4x4 (float3x3 m)
{
    return float4x4 (
        m[0].xyz, 0,
        m[1].xyz, 0,
        m[2].xyz, 0,
         0, 0, 0, 1
    );
}

#define I2x2 float2x2( \
    1, 0, \
    0, 1  \
)

#define I3x3 float3x3( \
    1, 0, 0, \
    0, 1, 0, \
    0, 0, 1  \
)

#define I4x4 float4x4( \
    1, 0, 0, 0, \
    0, 1, 0, 0, \
    0, 0, 1, 0, \
    0, 0, 0, 1  \
)

float2x2 diag (float2 v)
{
    return float2x2(
        v.x,   0,
          0, v.y
    );
}

float3x3 diag (float3 v)
{
    return float3x3(
        v.x,   0,   0,
          0, v.y,   0,
          0,   0, v.z
    );
}

float4x4 diag (float4 v)
{
    return float4x4(
        v.x,   0,   0,   0,
          0, v.y,   0,   0,
          0,   0, v.z,   0,
          0,   0,   0, v.w
    );
}

float3 zinv (float3 v)
{
    return float3(v.xy, -v.z);
}

float4 zinv (float4 v)
{
    return float4(v.xy, -v.z, v.w);
}

#define ZINV_MATRIX float3x3( \
    1,  0,  0, \
    0,  1,  0, \
    0,  0, -1  \
)

#define ZINV_MATRIX4 float4x4( \
    1,  0,  0,  0, \
    0,  1,  0,  0, \
    0,  0, -1,  0, \
    0,  0,  0,  1  \
)

float3x3 rotation_tri (float3 axis, float cos_, float sin_)
{
    float3x3 r = mul((float3x1)axis, (float1x3)axis) * (1 - cos_);
    r += float3x3(
        1, 0, 0,
        0, 1, 0,
        0, 0, 1
    ) * cos_;
    r += float3x3(
              0, -axis.z,  axis.y,
         axis.z,       0, -axis.x,
        -axis.y,  axis.x,       0
    ) * sin_;
    return r;
}

float3x3 rotation (float3 axis, float angle)
{
    return rotation_tri(axis, cos(angle), sin(angle));
}

float3x3 rotation (float3 theta)
{
    float angle = length(theta);
    return rotation_tri(theta / angle, cos(angle), sin(angle));
}

float4x4 rotation4_tri (float3 axis, float cos_, float sin_)
{
    float4 n = float4(axis, 0);
    float4x4 r = mul((float4x1)n, (float1x4)n) * (1 - cos_);
    r += float4x4(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 0
    ) * cos_;
    r += float4x4(
           0, -n.z,  n.y,    0,
         n.z,    0, -n.x,    0,
        -n.y,  n.x,    0,    0,
           0,    0,    0,    0
    ) * sin_;
    r[3][3] = 1;
    return r;
}

float4x4 rotaiton4 (float3 axis, float angle)
{
    return rotation4_tri(axis, cos(angle), sin(angle));
}

float4x4 rotation4 (float3 theta)
{
    float angle = length(theta);
    return rotation4_tri(theta / angle, cos(angle), sin(angle));
}

float clip2depth (float4 clip_pos)
{
    #if UNITY_UV_STARTS_AT_TOP
        return clip_pos.z / clip_pos.w;
    #else
        return clip_pos.z / clip_pos.w * 0.5 + 0.5;
    #endif
}

sampler3D _DitherMaskLOD;

float dither (float2 uv, float alpha)
{
    return tex3Dlod(_DitherMaskLOD, float4(frac(uv / 4.0), 15.0 / 16.0 * alpha, 0)).a;
}

float gaussian (float sigma, float x)
{
    x /= sigma;
    return exp(- 0.5 * x * x);
}

// sigma: (0.5, 0.8333...) if in 3 sigma
float3 gaussian_filter5 (float sigma)
{
    float3 filter;
    [unroll]
    for (int i = 0; i < 3; i++) {
        filter[i] = gaussian(sigma, i);
    }
    return filter / ((filter[2] + filter[1]) * 2 + filter[0]);
}

float3 xblur (sampler2D tex, float3 filter, float2 uv, float texel_size, float threshold=0.0, float supress=0.0)
{
    float3 color = 0;
    float3 tex_color;
    [unroll]
    for (int i = -2; i <= 2; i++) {
        tex_color = tex2D(tex, uv + i * float2(texel_size, 0)).xyz;
        color += any(tex_color > threshold) ? filter[abs(i)] * saturate(tex_color - supress) : 0;
    }
    return color;
}

float3 yblur (sampler2D tex, float3 filter, float2 uv, float texel_size, float threshold=0.0, float supress=0.0)
{
    float3 color = 0;
    float3 tex_color;
    [unroll]
    for (int i = -2; i <= 2; i++) {
        tex_color = tex2D(tex, uv + i * float2(0, texel_size)).xyz;
        color += any(tex_color > threshold) ? filter[abs(i)] * saturate(tex_color - supress) : 0;
    }
    return color;
}

#endif