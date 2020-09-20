% Calculate Gravitation load for parallel robot
% P3RRRRR10V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d1,d2,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 22:09
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR10V1G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:38:21
% EndTime: 2020-08-06 21:38:26
% DurationCPUTime: 5.19s
% Computational Cost: add. (1566->341), mult. (3708->570), div. (54->13), fcn. (2709->26), ass. (0->246)
t3094 = cos(pkin(3));
t3083 = t3094 ^ 2;
t3261 = pkin(6) * t3083;
t3098 = sin(qJ(3,3));
t3267 = pkin(2) * t3098;
t3059 = pkin(1) * t3267;
t3099 = sin(qJ(2,3));
t3070 = t3099 * pkin(6);
t3093 = sin(pkin(3));
t3107 = cos(qJ(3,3));
t3108 = cos(qJ(2,3));
t3223 = t3107 * t3108;
t3264 = pkin(2) * t3107;
t3202 = t3099 * t3264;
t3030 = -t3108 * pkin(6) + t3202;
t3245 = t3030 * t3094;
t2988 = 0.1e1 / (pkin(1) * t3245 + (-t3059 + (pkin(2) * t3223 + t3070) * pkin(5)) * t3093);
t3095 = legFrame(3,3);
t3062 = sin(t3095);
t3065 = cos(t3095);
t3024 = -g(1) * t3062 + g(2) * t3065;
t3027 = g(1) * t3065 + g(2) * t3062;
t3100 = sin(qJ(1,3));
t3109 = cos(qJ(1,3));
t3273 = m(3) * rSges(3,2);
t3274 = m(3) * rSges(3,1);
t3275 = m(2) * rSges(2,1);
t3021 = -t3098 * t3273 + t3107 * t3274 + t3275;
t3041 = m(2) * rSges(2,2) - rSges(3,3) * m(3);
t3197 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t3129 = t3021 * t3108 - t3041 * t3099 + t3197;
t3213 = m(2) * (rSges(2,3) + pkin(5)) + pkin(5) * m(3);
t3276 = m(1) * rSges(1,2);
t3135 = (t3021 * t3099 + t3041 * t3108) * t3094 - ((rSges(3,1) * t3098 + rSges(3,2) * t3107) * m(3) + t3213) * t3093 + t3276;
t3254 = ((t3129 * t3100 + t3135 * t3109) * t3027 + (t3135 * t3100 - t3129 * t3109) * t3024) * t2988;
t3101 = sin(qJ(3,2));
t3266 = pkin(2) * t3101;
t3060 = pkin(1) * t3266;
t3102 = sin(qJ(2,2));
t3071 = t3102 * pkin(6);
t3110 = cos(qJ(3,2));
t3111 = cos(qJ(2,2));
t3220 = t3110 * t3111;
t3263 = pkin(2) * t3110;
t3200 = t3102 * t3263;
t3031 = -t3111 * pkin(6) + t3200;
t3244 = t3031 * t3094;
t2989 = 0.1e1 / (pkin(1) * t3244 + (-t3060 + (pkin(2) * t3220 + t3071) * pkin(5)) * t3093);
t3096 = legFrame(2,3);
t3063 = sin(t3096);
t3066 = cos(t3096);
t3025 = -g(1) * t3063 + g(2) * t3066;
t3028 = g(1) * t3066 + g(2) * t3063;
t3103 = sin(qJ(1,2));
t3112 = cos(qJ(1,2));
t3022 = -t3101 * t3273 + t3110 * t3274 + t3275;
t3128 = t3022 * t3111 - t3041 * t3102 + t3197;
t3134 = (t3022 * t3102 + t3041 * t3111) * t3094 - ((rSges(3,1) * t3101 + rSges(3,2) * t3110) * m(3) + t3213) * t3093 + t3276;
t3253 = ((t3128 * t3103 + t3134 * t3112) * t3028 + (t3134 * t3103 - t3128 * t3112) * t3025) * t2989;
t3104 = sin(qJ(3,1));
t3265 = pkin(2) * t3104;
t3061 = pkin(1) * t3265;
t3105 = sin(qJ(2,1));
t3072 = t3105 * pkin(6);
t3113 = cos(qJ(3,1));
t3114 = cos(qJ(2,1));
t3217 = t3113 * t3114;
t3262 = pkin(2) * t3113;
t3198 = t3105 * t3262;
t3032 = -t3114 * pkin(6) + t3198;
t3243 = t3032 * t3094;
t2990 = 0.1e1 / (pkin(1) * t3243 + (-t3061 + (pkin(2) * t3217 + t3072) * pkin(5)) * t3093);
t3097 = legFrame(1,3);
t3064 = sin(t3097);
t3067 = cos(t3097);
t3026 = -g(1) * t3064 + g(2) * t3067;
t3029 = g(1) * t3067 + g(2) * t3064;
t3106 = sin(qJ(1,1));
t3115 = cos(qJ(1,1));
t3023 = -t3104 * t3273 + t3113 * t3274 + t3275;
t3127 = t3023 * t3114 - t3041 * t3105 + t3197;
t3133 = (t3023 * t3105 + t3041 * t3114) * t3094 - ((rSges(3,1) * t3104 + rSges(3,2) * t3113) * m(3) + t3213) * t3093 + t3276;
t3252 = ((t3127 * t3106 + t3133 * t3115) * t3029 + (t3133 * t3106 - t3127 * t3115) * t3026) * t2990;
t3236 = t3093 * t3094;
t3148 = t3026 * t3115 - t3029 * t3106;
t3260 = g(3) * t3093;
t3303 = t3148 * t3094 + t3260;
t3149 = t3025 * t3112 - t3028 * t3103;
t3302 = t3149 * t3094 + t3260;
t3150 = t3024 * t3109 - t3027 * t3100;
t3301 = t3150 * t3094 + t3260;
t3224 = t3104 * t3114;
t3172 = t3105 * t3224;
t3210 = t3093 * t3262;
t3090 = t3113 ^ 2;
t3215 = t3114 * t3090;
t3233 = t3093 * t3104;
t3290 = t3072 + pkin(1);
t3300 = t3094 * (pkin(2) * t3215 + t3113 * t3290) - t3172 * t3210 + pkin(6) * (t3114 - 0.1e1) * (t3114 + 0.1e1) * t3233;
t3225 = t3101 * t3111;
t3173 = t3102 * t3225;
t3211 = t3093 * t3263;
t3087 = t3110 ^ 2;
t3218 = t3111 * t3087;
t3234 = t3093 * t3101;
t3291 = t3071 + pkin(1);
t3299 = t3094 * (pkin(2) * t3218 + t3110 * t3291) - t3173 * t3211 + pkin(6) * (t3111 - 0.1e1) * (t3111 + 0.1e1) * t3234;
t3226 = t3098 * t3108;
t3174 = t3099 * t3226;
t3212 = t3093 * t3264;
t3084 = t3107 ^ 2;
t3221 = t3108 * t3084;
t3235 = t3093 * t3098;
t3292 = t3070 + pkin(1);
t3298 = t3094 * (pkin(2) * t3221 + t3107 * t3292) - t3174 * t3212 + pkin(6) * (t3108 - 0.1e1) * (t3108 + 0.1e1) * t3235;
t3297 = -0.2e1 * pkin(6);
t3296 = 0.2e1 * pkin(6);
t3295 = pkin(1) * t3099;
t3294 = pkin(1) * t3102;
t3293 = pkin(1) * t3105;
t3259 = g(3) * t3094;
t3017 = t3064 * t3106 - t3067 * t3115;
t3118 = pkin(6) ^ 2;
t3119 = pkin(2) ^ 2;
t3166 = t3090 * t3119 - t3118;
t3289 = t3166 * t3017;
t3016 = t3063 * t3103 - t3066 * t3112;
t3167 = t3087 * t3119 - t3118;
t3288 = t3167 * t3016;
t3015 = t3062 * t3100 - t3065 * t3109;
t3168 = t3084 * t3119 - t3118;
t3287 = t3168 * t3015;
t3086 = t3108 ^ 2;
t3036 = (t3086 - 0.2e1) * t3267 - pkin(5);
t3043 = pkin(5) * t3098 + pkin(2);
t3270 = pkin(2) * t3084;
t3171 = t3043 - 0.2e1 * t3270;
t3280 = (-pkin(6) * t3174 - t3036 * t3107) * t3236 + t3223 * t3261 + (t3171 * t3083 - t3043 + t3270) * t3099;
t3089 = t3111 ^ 2;
t3037 = (t3089 - 0.2e1) * t3266 - pkin(5);
t3046 = pkin(5) * t3101 + pkin(2);
t3269 = pkin(2) * t3087;
t3170 = t3046 - 0.2e1 * t3269;
t3279 = (-pkin(6) * t3173 - t3037 * t3110) * t3236 + t3220 * t3261 + (t3170 * t3083 - t3046 + t3269) * t3102;
t3092 = t3114 ^ 2;
t3038 = (t3092 - 0.2e1) * t3265 - pkin(5);
t3049 = pkin(5) * t3104 + pkin(2);
t3268 = pkin(2) * t3090;
t3169 = t3049 - 0.2e1 * t3268;
t3278 = (-pkin(6) * t3172 - t3038 * t3113) * t3236 + t3217 * t3261 + (t3169 * t3083 - t3049 + t3268) * t3105;
t3277 = (t3083 - 0.1e1) * t3296;
t3272 = m(3) / pkin(2);
t3271 = pkin(1) * t3094;
t3255 = rSges(3,2) * t3093;
t2991 = t3024 * t3100 + t3027 * t3109;
t3039 = t3041 * t3260;
t3085 = 0.1e1 / t3107;
t3222 = t3108 * t2991;
t3229 = t3094 * t3099;
t3251 = (t3039 * t3099 + (t3150 * t3229 + t3222) * t3041 + (t2991 * t3099 - t3301 * t3108) * t3021) * t3085;
t2992 = t3025 * t3103 + t3028 * t3112;
t3088 = 0.1e1 / t3110;
t3219 = t3111 * t2992;
t3228 = t3094 * t3102;
t3250 = (t3039 * t3102 + (t3149 * t3228 + t3219) * t3041 + (t2992 * t3102 - t3302 * t3111) * t3022) * t3088;
t2993 = t3026 * t3106 + t3029 * t3115;
t3091 = 0.1e1 / t3113;
t3216 = t3114 * t2993;
t3227 = t3094 * t3105;
t3249 = (t3039 * t3105 + (t3148 * t3227 + t3216) * t3041 + (t2993 * t3105 - t3303 * t3114) * t3023) * t3091;
t3044 = pkin(5) + t3267;
t3242 = t3044 * t3093;
t3047 = pkin(5) + t3266;
t3241 = t3047 * t3093;
t3050 = pkin(5) + t3265;
t3240 = t3050 * t3093;
t3239 = t3093 * (-pkin(5) * t3070 + t3059);
t3238 = t3093 * (-pkin(5) * t3071 + t3060);
t3237 = t3093 * (-pkin(5) * t3072 + t3061);
t3232 = t3093 * t3108;
t3231 = t3093 * t3111;
t3230 = t3093 * t3114;
t3214 = rSges(3,2) * t3259;
t3055 = pkin(6) * t3271;
t3209 = t3094 * t3264;
t3208 = t3094 * t3263;
t3207 = t3094 * t3262;
t3052 = pkin(6) + t3295;
t3053 = pkin(6) + t3294;
t3054 = pkin(6) + t3293;
t3193 = t2988 * t3251;
t3192 = t2989 * t3250;
t3191 = t2990 * t3249;
t3190 = t3015 * t3235;
t3189 = t3016 * t3234;
t3188 = t3017 * t3233;
t3018 = t3062 * t3109 + t3065 * t3100;
t3187 = t3018 * t3235;
t3019 = t3063 * t3112 + t3066 * t3103;
t3186 = t3019 * t3234;
t3020 = t3064 * t3115 + t3067 * t3106;
t3185 = t3020 * t3233;
t3184 = t3100 * t3242;
t3183 = t3109 * t3242;
t3182 = t3103 * t3241;
t3181 = t3112 * t3241;
t3180 = t3106 * t3240;
t3179 = t3115 * t3240;
t3178 = (t3094 + 0.1e1) * (t3094 - 0.1e1) * t3119;
t3177 = t3094 * t3232;
t3176 = t3094 * t3231;
t3175 = t3094 * t3230;
t3123 = t3301 * t3099 + t3222;
t3165 = (((t3150 * t3093 - t3259) * rSges(3,1) + t3123 * rSges(3,2)) * t3107 + t3098 * (t3123 * rSges(3,1) - t3150 * t3255 + t3214)) * t3085 * t3272;
t3122 = t3302 * t3102 + t3219;
t3164 = (((t3149 * t3093 - t3259) * rSges(3,1) + t3122 * rSges(3,2)) * t3110 + t3101 * (t3122 * rSges(3,1) - t3149 * t3255 + t3214)) * t3088 * t3272;
t3121 = t3303 * t3105 + t3216;
t3163 = (((t3148 * t3093 - t3259) * rSges(3,1) + t3121 * rSges(3,2)) * t3113 + t3104 * (t3121 * rSges(3,1) - t3148 * t3255 + t3214)) * t3091 * t3272;
t3162 = t3015 * t3209;
t3161 = t3016 * t3208;
t3160 = t3017 * t3207;
t3159 = t3018 * t3209;
t3158 = t3019 * t3208;
t3157 = t3020 * t3207;
t3153 = t3168 * t3018;
t3152 = t3167 * t3019;
t3151 = t3166 * t3020;
t3147 = 0.1e1 / ((-pkin(5) * t3212 + t3055) * t3108 - t3202 * t3271 + t3239) * t3093 * t3165;
t3146 = 0.1e1 / ((-pkin(5) * t3211 + t3055) * t3111 - t3200 * t3271 + t3238) * t3093 * t3164;
t3145 = 0.1e1 / ((-pkin(5) * t3210 + t3055) * t3114 - t3198 * t3271 + t3237) * t3093 * t3163;
t3045 = 0.2e1 * t3070 + pkin(1);
t3138 = t3045 * t3109 + t3184;
t3048 = 0.2e1 * t3071 + pkin(1);
t3137 = t3048 * t3112 + t3182;
t3051 = 0.2e1 * t3072 + pkin(1);
t3136 = t3051 * t3115 + t3180;
t3132 = t3052 * t3109 + t3099 * t3184;
t3131 = t3053 * t3112 + t3102 * t3182;
t3130 = t3054 * t3115 + t3105 * t3180;
t3008 = t3051 * t3106 - t3179;
t3007 = t3048 * t3103 - t3181;
t3006 = t3045 * t3100 - t3183;
t3002 = t3054 * t3106 - t3105 * t3179;
t3001 = t3053 * t3103 - t3102 * t3181;
t3000 = t3052 * t3100 - t3099 * t3183;
t1 = [-(t3017 * t3072 + t3020 * t3243 + (t3017 * t3217 - t3185) * pkin(2)) * t3252 - (-t3300 * t3017 + t3278 * t3020 + t3188 * t3293) * t3191 + ((t3157 * t3297 + t3289) * t3092 + ((t3064 * t3008 - t3136 * t3067) * t3262 + t3151 * t3227) * t3114 + (t3064 * t3002 - t3130 * t3067 + t3157) * pkin(6)) * t3145 - (t3016 * t3071 + t3019 * t3244 + (t3016 * t3220 - t3186) * pkin(2)) * t3253 - (-t3299 * t3016 + t3279 * t3019 + t3189 * t3294) * t3192 + ((t3158 * t3297 + t3288) * t3089 + ((t3063 * t3007 - t3137 * t3066) * t3263 + t3152 * t3228) * t3111 + (t3063 * t3001 - t3131 * t3066 + t3158) * pkin(6)) * t3146 - (t3015 * t3070 + t3018 * t3245 + (t3015 * t3223 - t3187) * pkin(2)) * t3254 - (-t3298 * t3015 + t3280 * t3018 + t3190 * t3295) * t3193 + ((t3159 * t3297 + t3287) * t3086 + ((t3062 * t3006 - t3138 * t3065) * t3264 + t3153 * t3229) * t3108 + (t3062 * t3000 - t3132 * t3065 + t3159) * pkin(6)) * t3147 - m(4) * g(1); -(-t3020 * t3072 + t3017 * t3243 + (-t3020 * t3217 - t3188) * pkin(2)) * t3252 - (t3278 * t3017 + t3300 * t3020 - t3185 * t3293) * t3191 - ((t3160 * t3296 + t3151) * t3092 + ((t3008 * t3067 + t3064 * t3136) * t3262 - t3227 * t3289) * t3114 + (t3002 * t3067 + t3064 * t3130 - t3160) * pkin(6)) * t3145 - (-t3019 * t3071 + t3016 * t3244 + (-t3019 * t3220 - t3189) * pkin(2)) * t3253 - (t3279 * t3016 + t3299 * t3019 - t3186 * t3294) * t3192 - ((t3161 * t3296 + t3152) * t3089 + ((t3007 * t3066 + t3063 * t3137) * t3263 - t3228 * t3288) * t3111 + (t3001 * t3066 + t3063 * t3131 - t3161) * pkin(6)) * t3146 - (-t3018 * t3070 + t3015 * t3245 + (-t3018 * t3223 - t3190) * pkin(2)) * t3254 - (t3280 * t3015 + t3298 * t3018 - t3187 * t3295) * t3193 - ((t3162 * t3296 + t3153) * t3086 + ((t3006 * t3065 + t3062 * t3138) * t3264 - t3229 * t3287) * t3108 + (t3000 * t3065 + t3062 * t3132 - t3162) * pkin(6)) * t3147 - m(4) * g(2); (t3032 * t3093 + t3094 * t3265) * t3252 + (t3031 * t3093 + t3094 * t3266) * t3253 + (t3030 * t3093 + t3094 * t3267) * t3254 - m(4) * g(3) + (((pkin(6) * t3175 + t3038 * t3083 + pkin(5) + (-t3092 + 0.1e1) * t3265) * t3113 - t3290 * t3224 + (t3169 * t3236 + t3224 * t3261) * t3105) * t3249 + (-t3105 * t3178 * t3215 + (t3050 * t3175 + t3092 * t3277 + t3054 - t3261) * t3262 - ((-t3083 * t3072 + t3290) * t3114 - t3227 * t3240) * pkin(6)) * t3163) / ((pkin(1) * t3227 + pkin(5) * t3230) * t3262 - t3114 * t3055 - t3237) + (((pkin(6) * t3176 + t3037 * t3083 + pkin(5) + (-t3089 + 0.1e1) * t3266) * t3110 - t3291 * t3225 + (t3170 * t3236 + t3225 * t3261) * t3102) * t3250 + (-t3102 * t3178 * t3218 + (t3047 * t3176 + t3089 * t3277 + t3053 - t3261) * t3263 - ((-t3083 * t3071 + t3291) * t3111 - t3228 * t3241) * pkin(6)) * t3164) / ((pkin(1) * t3228 + pkin(5) * t3231) * t3263 - t3111 * t3055 - t3238) + (((pkin(6) * t3177 + t3036 * t3083 + pkin(5) + (-t3086 + 0.1e1) * t3267) * t3107 - t3292 * t3226 + (t3171 * t3236 + t3226 * t3261) * t3099) * t3251 + (-t3099 * t3178 * t3221 + (t3044 * t3177 + t3086 * t3277 + t3052 - t3261) * t3264 - ((-t3083 * t3070 + t3292) * t3108 - t3229 * t3242) * pkin(6)) * t3165) / ((pkin(1) * t3229 + pkin(5) * t3232) * t3264 - t3108 * t3055 - t3239);];
taugX  = t1;
