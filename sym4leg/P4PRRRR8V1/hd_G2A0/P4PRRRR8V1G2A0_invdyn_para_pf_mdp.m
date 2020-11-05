% Calculate minimal parameter regressor of inverse dynamics forces for
% P4PRRRR8V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% xDDP [4x1]
%   Generalized platform accelerations
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P4PRRRR8V1G2A0_convert_par2_MPV_fixb.m

% Output:
% tauX [4x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:06
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR8V1G2A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(6,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_mdp: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_mdp: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:05:00
% EndTime: 2020-08-07 11:05:23
% DurationCPUTime: 23.78s
% Computational Cost: add. (119159->717), mult. (280638->1344), div. (11524->18), fcn. (219964->30), ass. (0->543)
t3081 = cos(qJ(3,4));
t3066 = 0.1e1 / t3081;
t3082 = cos(qJ(2,4));
t3080 = sin(qJ(2,4));
t3306 = t3080 * t3081;
t3029 = pkin(2) * t3306 - pkin(5) * t3082;
t3076 = sin(pkin(3));
t3078 = cos(pkin(3));
t3079 = sin(qJ(3,4));
t3405 = pkin(2) * t3079;
t3232 = t3029 * t3076 + t3078 * t3405;
t3406 = 0.1e1 / t3232;
t3349 = t3406 * t3066;
t3101 = cos(qJ(3,1));
t3072 = 0.1e1 / t3101;
t3102 = cos(qJ(2,1));
t3096 = sin(qJ(2,1));
t3297 = t3096 * t3101;
t3033 = pkin(2) * t3297 - pkin(5) * t3102;
t3095 = sin(qJ(3,1));
t3401 = pkin(2) * t3095;
t3229 = t3033 * t3076 + t3078 * t3401;
t3409 = 0.1e1 / t3229;
t3340 = t3409 * t3072;
t3083 = legFrame(4,2);
t3056 = sin(t3083);
t3257 = t3056 * t3349;
t3086 = legFrame(1,2);
t3059 = sin(t3086);
t3242 = t3059 * t3340;
t3077 = cos(pkin(6));
t3055 = t3077 * g(3);
t3075 = sin(pkin(6));
t3060 = cos(t3083);
t3180 = g(1) * t3060 - g(2) * t3056;
t2997 = t3075 * t3180 + t3055;
t3107 = xP(4);
t3064 = sin(t3107);
t3065 = cos(t3107);
t3108 = koppelP(4,2);
t3112 = koppelP(4,1);
t3019 = t3064 * t3112 + t3065 * t3108;
t3023 = -t3064 * t3108 + t3065 * t3112;
t3103 = xDP(4);
t3105 = xDP(2);
t3106 = xDP(1);
t2917 = (-t3019 * t3103 + t3106) * t3060 - (t3023 * t3103 + t3105) * t3056;
t3104 = xDP(3);
t3325 = t3075 * t3104;
t2909 = t2917 * t3077 - t3325;
t3045 = t3104 * t3077;
t2913 = t2917 * t3075 + t3045;
t3318 = t3078 * t3082;
t3319 = t3078 * t3080;
t3404 = pkin(2) * t3081;
t2880 = -(t2909 * t3080 + t2913 * t3318) * t3404 + pkin(5) * (t2909 * t3082 - t2913 * t3319);
t3116 = pkin(5) ^ 2;
t3117 = pkin(2) ^ 2;
t3327 = t3066 * t3079;
t3255 = t3406 * t3327;
t3311 = t3078 * t3104;
t3324 = t3076 * t3081;
t2892 = ((-t2917 * t3319 - t3082 * t3104) * t3075 + (t2917 * t3082 - t3080 * t3311) * t3077) * t3079 - t2913 * t3324;
t3372 = t2892 * t3406;
t2863 = -pkin(5) * t2880 * t3255 + (t3066 * t3116 + t3081 * t3117) * t3372;
t3296 = pkin(5) * t3372;
t3228 = t3079 * t3296;
t3381 = t2880 * t3406;
t2867 = (t3228 - t3381) * t3066;
t3030 = pkin(5) * t3080 + t3082 * t3404;
t3139 = -t3029 * t3078 + t3076 * t3405;
t2919 = -t3030 * t3075 + t3077 * t3139;
t2899 = t2919 * t3056 + t3060 * t3232;
t2900 = -t2919 * t3060 + t3056 * t3232;
t2920 = t3030 * t3077 + t3075 * t3139;
t3074 = t3103 ^ 2;
t3087 = xDDP(4);
t3089 = xDDP(2);
t2947 = -t3019 * t3074 + t3023 * t3087 + t3089;
t3090 = xDDP(1);
t2951 = -t3019 * t3087 - t3023 * t3074 + t3090;
t3088 = xDDP(3);
t3281 = t3406 * t3372;
t2967 = 0.1e1 / t3232 ^ 2;
t3380 = t2880 * t2967;
t2848 = (t2899 * t2947 + t2900 * t2951 + t2920 * t3088) * t3406 + (t2863 * t3281 - t2867 * t3380) * t3066;
t3320 = t3077 * t3078;
t3027 = -t3076 * g(1) + g(2) * t3320;
t3028 = g(1) * t3320 + g(2) * t3076;
t3397 = g(3) * t3075;
t3040 = t3078 * t3397;
t3127 = t2848 * t3076 + t3027 * t3056 - t3028 * t3060 + t3040;
t2831 = t3080 * t2997 + t3082 * t3127;
t2977 = 0.1e1 / t3229 ^ 2;
t3100 = cos(qJ(2,2));
t3094 = sin(qJ(2,2));
t3099 = cos(qJ(3,2));
t3300 = t3094 * t3099;
t3032 = pkin(2) * t3300 - pkin(5) * t3100;
t3093 = sin(qJ(3,2));
t3402 = pkin(2) * t3093;
t3230 = t3032 * t3076 + t3078 * t3402;
t2975 = 0.1e1 / t3230 ^ 2;
t3408 = 0.1e1 / t3230;
t3098 = cos(qJ(2,3));
t3092 = sin(qJ(2,3));
t3097 = cos(qJ(3,3));
t3303 = t3092 * t3097;
t3031 = pkin(2) * t3303 - pkin(5) * t3098;
t3091 = sin(qJ(3,3));
t3403 = pkin(2) * t3091;
t3231 = t3031 * t3076 + t3078 * t3403;
t2973 = 0.1e1 / t3231 ^ 2;
t3407 = 0.1e1 / t3231;
t3123 = t3101 ^ 2;
t3073 = 0.1e1 / t3123;
t3111 = koppelP(1,2);
t3115 = koppelP(1,1);
t3022 = t3064 * t3115 + t3065 * t3111;
t3026 = -t3064 * t3111 + t3065 * t3115;
t3063 = cos(t3086);
t2916 = (-t3022 * t3103 + t3106) * t3063 - (t3026 * t3103 + t3105) * t3059;
t2912 = t2916 * t3075 + t3045;
t3315 = t3078 * t3096;
t3321 = t3076 * t3101;
t2898 = ((-t2916 * t3315 - t3102 * t3104) * t3075 + (t2916 * t3102 - t3096 * t3311) * t3077) * t3095 - t2912 * t3321;
t3369 = t2898 ^ 2 * t2977;
t3278 = t3073 * t3369;
t3339 = t3409 * t3095;
t3158 = t3278 * t3339;
t3014 = t3075 * t3315 - t3077 * t3102;
t3312 = t3078 * t3102;
t3398 = pkin(2) * t3101;
t3233 = pkin(5) * t3014 + (t3075 * t3312 + t3077 * t3096) * t3398;
t3132 = t3233 * t3158;
t3122 = t3099 ^ 2;
t3071 = 0.1e1 / t3122;
t3110 = koppelP(2,2);
t3114 = koppelP(2,1);
t3021 = t3064 * t3114 + t3065 * t3110;
t3025 = -t3064 * t3110 + t3065 * t3114;
t3085 = legFrame(2,2);
t3058 = sin(t3085);
t3062 = cos(t3085);
t2915 = (-t3021 * t3103 + t3106) * t3062 - (t3025 * t3103 + t3105) * t3058;
t2911 = t2915 * t3075 + t3045;
t3316 = t3078 * t3094;
t3322 = t3076 * t3099;
t2897 = ((-t2915 * t3316 - t3100 * t3104) * t3075 + (t2915 * t3100 - t3094 * t3311) * t3077) * t3093 - t2911 * t3322;
t3370 = t2897 ^ 2 * t2975;
t3279 = t3071 * t3370;
t3342 = t3408 * t3093;
t3159 = t3279 * t3342;
t3313 = t3078 * t3100;
t3399 = pkin(2) * t3099;
t3234 = -pkin(5) * (-t3075 * t3316 + t3077 * t3100) + (t3075 * t3313 + t3077 * t3094) * t3399;
t3133 = t3234 * t3159;
t3121 = t3097 ^ 2;
t3069 = 0.1e1 / t3121;
t3109 = koppelP(3,2);
t3113 = koppelP(3,1);
t3020 = t3064 * t3113 + t3065 * t3109;
t3024 = -t3064 * t3109 + t3065 * t3113;
t3084 = legFrame(3,2);
t3057 = sin(t3084);
t3061 = cos(t3084);
t2918 = (-t3020 * t3103 + t3106) * t3061 - (t3024 * t3103 + t3105) * t3057;
t2914 = t2918 * t3075 + t3045;
t3317 = t3078 * t3092;
t3323 = t3076 * t3097;
t2896 = ((-t2918 * t3317 - t3098 * t3104) * t3075 + (t2918 * t3098 - t3092 * t3311) * t3077) * t3091 - t2914 * t3323;
t3371 = t2896 ^ 2 * t2973;
t3280 = t3069 * t3371;
t3345 = t3407 * t3091;
t3160 = t3280 * t3345;
t3314 = t3078 * t3098;
t3400 = pkin(2) * t3097;
t3235 = -pkin(5) * (-t3075 * t3317 + t3077 * t3098) + (t3075 * t3314 + t3077 * t3092) * t3400;
t3134 = t3235 * t3160;
t3120 = t3081 ^ 2;
t3067 = 0.1e1 / t3120;
t3373 = t2892 ^ 2 * t2967;
t3282 = t3067 * t3373;
t3348 = t3406 * t3079;
t3161 = t3282 * t3348;
t3013 = -t3075 * t3319 + t3077 * t3082;
t3236 = -pkin(5) * t3013 + (t3075 * t3318 + t3077 * t3080) * t3404;
t3135 = t3236 * t3161;
t3396 = t3089 - g(2);
t3395 = t3090 - g(1);
t3118 = 0.1e1 / pkin(2);
t3394 = MDP(9) * t3118;
t3309 = t3079 * t3080;
t3131 = t3078 * t3309 + t3324;
t3308 = t3079 * t3082;
t2955 = t3075 * t3131 - t3077 * t3308;
t3347 = t3406 * t3118;
t3254 = t3066 * t3347;
t3168 = t2880 * t3076 * t3254;
t3256 = t3060 * t3349;
t3196 = t2955 * t3256;
t2956 = -t3075 * t3308 - t3077 * t3131;
t3262 = t2956 * t3349;
t3310 = t3078 * t3118;
t2809 = t2955 * t2947 * t3257 - t2951 * t3196 + t3088 * t3262 + (-(t3078 * t2867 + (pkin(2) * (t3076 * t3082 * t3372 + t3310 * t3381) * t3120 - t3076 * (t2880 * t3348 - t3296) * t3306) * t3066) * t3281 - (t3082 * t3168 + (-t3076 * t3309 + (-t3066 + t3081) * t3078) * t3372) * t3380) * t3067;
t3393 = t2809 * t3406;
t3392 = t2809 * t3079;
t3391 = t2809 * t3082;
t3068 = 0.1e1 / t3097;
t3368 = t2896 * t3407;
t3295 = pkin(5) * t3368;
t3227 = t3091 * t3295;
t2910 = t2918 * t3077 - t3325;
t2884 = -(t2910 * t3092 + t2914 * t3314) * t3400 + pkin(5) * (t2910 * t3098 - t2914 * t3317);
t3379 = t2884 * t3407;
t2869 = (t3227 - t3379) * t3068;
t2948 = -t3020 * t3074 + t3024 * t3087 + t3089;
t2952 = -t3020 * t3087 - t3024 * t3074 + t3090;
t3344 = t3407 * t3118;
t3249 = t3068 * t3344;
t3166 = t2884 * t3076 * t3249;
t3305 = t3091 * t3092;
t3130 = t3078 * t3305 + t3323;
t3304 = t3091 * t3098;
t2957 = t3075 * t3130 - t3077 * t3304;
t3346 = t3407 * t3068;
t3251 = t3061 * t3346;
t3193 = t2957 * t3251;
t3252 = t3057 * t3346;
t3194 = t2957 * t3252;
t2960 = -t3075 * t3304 - t3077 * t3130;
t3261 = t2960 * t3346;
t3277 = t3407 * t3368;
t3378 = t2884 * t2973;
t2810 = t2948 * t3194 - t2952 * t3193 + t3088 * t3261 + (-(t3078 * t2869 + (pkin(2) * (t3076 * t3098 * t3368 + t3310 * t3379) * t3121 - t3076 * (t2884 * t3345 - t3295) * t3303) * t3068) * t3277 + (-t3098 * t3166 + (t3076 * t3305 + (t3068 - t3097) * t3078) * t3368) * t3378) * t3069;
t3390 = t2810 * t3407;
t3389 = t2810 * t3091;
t3388 = t2810 * t3098;
t3070 = 0.1e1 / t3099;
t3367 = t2897 * t3408;
t3294 = pkin(5) * t3367;
t3226 = t3093 * t3294;
t2907 = t2915 * t3077 - t3325;
t2885 = -(t2907 * t3094 + t2911 * t3313) * t3399 + pkin(5) * (t2907 * t3100 - t2911 * t3316);
t3377 = t2885 * t3408;
t2870 = (t3226 - t3377) * t3070;
t2949 = -t3021 * t3074 + t3025 * t3087 + t3089;
t2953 = -t3021 * t3087 - t3025 * t3074 + t3090;
t3341 = t3408 * t3118;
t3244 = t3070 * t3341;
t3164 = t2885 * t3076 * t3244;
t3302 = t3093 * t3094;
t3129 = t3078 * t3302 + t3322;
t3301 = t3093 * t3100;
t2958 = t3075 * t3129 - t3077 * t3301;
t3343 = t3408 * t3070;
t3246 = t3062 * t3343;
t3190 = t2958 * t3246;
t3247 = t3058 * t3343;
t3191 = t2958 * t3247;
t2961 = -t3075 * t3301 - t3077 * t3129;
t3260 = t2961 * t3343;
t3276 = t3408 * t3367;
t3376 = t2885 * t2975;
t2811 = t2949 * t3191 - t2953 * t3190 + t3088 * t3260 + (-(t3078 * t2870 + (pkin(2) * (t3076 * t3100 * t3367 + t3310 * t3377) * t3122 - t3076 * (t2885 * t3342 - t3294) * t3300) * t3070) * t3276 + (-t3100 * t3164 + (t3076 * t3302 + (t3070 - t3099) * t3078) * t3367) * t3376) * t3071;
t3387 = t2811 * t3408;
t3386 = t2811 * t3093;
t3385 = t2811 * t3100;
t3366 = t2898 * t3409;
t3293 = pkin(5) * t3366;
t3225 = t3095 * t3293;
t2908 = t2916 * t3077 - t3325;
t2886 = -(t2908 * t3096 + t2912 * t3312) * t3398 + pkin(5) * (t2908 * t3102 - t2912 * t3315);
t3375 = t2886 * t3409;
t2871 = (t3225 - t3375) * t3072;
t2950 = -t3022 * t3074 + t3026 * t3087 + t3089;
t2954 = -t3022 * t3087 - t3026 * t3074 + t3090;
t3299 = t3095 * t3096;
t3128 = t3078 * t3299 + t3321;
t3298 = t3095 * t3102;
t2959 = t3075 * t3128 - t3077 * t3298;
t3338 = t3409 * t3118;
t3239 = t3072 * t3338;
t3162 = t2886 * t3076 * t3239;
t3241 = t3063 * t3340;
t3188 = t2959 * t3241;
t2962 = -t3075 * t3298 - t3077 * t3128;
t3259 = t2962 * t3340;
t3275 = t3409 * t3366;
t3374 = t2886 * t2977;
t2812 = t2959 * t2950 * t3242 - t2954 * t3188 + t3088 * t3259 + (-(t3078 * t2871 + (pkin(2) * (t3076 * t3102 * t3366 + t3310 * t3375) * t3123 - t3076 * (t2886 * t3339 - t3293) * t3297) * t3072) * t3275 + (-t3102 * t3162 + (t3076 * t3299 + (t3072 - t3101) * t3078) * t3366) * t3374) * t3073;
t3384 = t2812 * t3409;
t3383 = t2812 * t3095;
t3382 = t2812 * t3102;
t3365 = t2899 * t3406;
t3364 = t2900 * t3406;
t3034 = pkin(5) * t3092 + t3098 * t3400;
t3138 = -t3031 * t3078 + t3076 * t3403;
t2921 = -t3034 * t3075 + t3077 * t3138;
t2901 = t2921 * t3057 + t3061 * t3231;
t3363 = t2901 * t3407;
t2902 = -t2921 * t3061 + t3057 * t3231;
t3362 = t2902 * t3407;
t3035 = pkin(5) * t3094 + t3100 * t3399;
t3137 = -t3032 * t3078 + t3076 * t3402;
t2922 = -t3035 * t3075 + t3077 * t3137;
t2903 = t2922 * t3058 + t3062 * t3230;
t3361 = t2903 * t3408;
t2904 = -t2922 * t3062 + t3058 * t3230;
t3360 = t2904 * t3408;
t3036 = pkin(5) * t3096 + t3102 * t3398;
t3136 = -t3033 * t3078 + t3076 * t3401;
t2923 = -t3036 * t3075 + t3077 * t3136;
t2905 = t2923 * t3059 + t3063 * t3229;
t3359 = t2905 * t3409;
t2906 = -t2923 * t3063 + t3059 * t3229;
t3358 = t2906 * t3409;
t3357 = t2920 * t3406;
t2924 = t3034 * t3077 + t3075 * t3138;
t3356 = t2924 * t3407;
t2925 = t3035 * t3077 + t3075 * t3137;
t3355 = t2925 * t3408;
t2926 = t3036 * t3077 + t3075 * t3136;
t3354 = t2926 * t3409;
t3353 = t2955 * t3406;
t3352 = t2957 * t3407;
t3351 = t2958 * t3408;
t3350 = t2959 * t3409;
t2931 = -(-t3075 * t3080 + t3077 * t3318) * t3404 - pkin(5) * (t3075 * t3082 + t3077 * t3319);
t3258 = t3406 * t3347;
t3216 = t2892 * t3258;
t2825 = (-t3078 * t2863 * t3216 - (-t3079 * t3029 * t3168 + t3078 * (-t3066 * t3228 + t3081 * t3381)) * t2880 * t3258) * t3067 + (t2931 * t3088 + (t2947 * t3056 - t2951 * t3060) * t3236) * t3254;
t3169 = t2880 * t3216;
t3143 = 0.2e1 * t3169;
t2813 = t3067 * t3082 * t3143 + t2825 * t3080;
t3119 = 0.1e1 / pkin(2) ^ 2;
t3286 = t2880 ^ 2 * t2967 * t3119;
t2817 = -t3067 * t3079 * t3286 + t2825 * t3081;
t3151 = -(t3286 + t3373) * t3067 * t3080 + t3391;
t2777 = (-t3079 * t2813 + t3081 * t3151) * t3076 + t3078 * t2817;
t3337 = t3406 * t2777;
t2818 = t2825 * t3079 + t3066 * t3286;
t2778 = (-t3081 * t2813 - t3079 * t3151) * t3076 - t3078 * t2818;
t3336 = t3406 * t2778;
t3250 = t3068 * t3345;
t2864 = -pkin(5) * t2884 * t3250 + (t3068 * t3116 + t3097 * t3117) * t3368;
t2934 = -(-t3075 * t3092 + t3077 * t3314) * t3400 - pkin(5) * (t3075 * t3098 + t3077 * t3317);
t3253 = t3407 * t3344;
t3215 = t2896 * t3253;
t2826 = (-t3078 * t2864 * t3215 - (-t3091 * t3031 * t3166 + t3078 * (-t3068 * t3227 + t3097 * t3379)) * t2884 * t3253) * t3069 + (t2934 * t3088 + (t2948 * t3057 - t2952 * t3061) * t3235) * t3249;
t3167 = t2884 * t3215;
t3142 = 0.2e1 * t3167;
t2814 = t3069 * t3098 * t3142 + t2826 * t3092;
t3285 = t2884 ^ 2 * t2973 * t3119;
t2819 = -t3069 * t3091 * t3285 + t2826 * t3097;
t3150 = -(t3285 + t3371) * t3069 * t3092 + t3388;
t2779 = (-t3091 * t2814 + t3097 * t3150) * t3076 + t3078 * t2819;
t3335 = t3407 * t2779;
t2822 = t2826 * t3091 + t3068 * t3285;
t2782 = (-t3097 * t2814 - t3091 * t3150) * t3076 - t3078 * t2822;
t3334 = t3407 * t2782;
t3245 = t3070 * t3342;
t2865 = -pkin(5) * t2885 * t3245 + (t3070 * t3116 + t3099 * t3117) * t3367;
t2935 = -(-t3075 * t3094 + t3077 * t3313) * t3399 - pkin(5) * (t3075 * t3100 + t3077 * t3316);
t3248 = t3408 * t3341;
t3214 = t2897 * t3248;
t2827 = (-t3078 * t2865 * t3214 - (-t3093 * t3032 * t3164 + t3078 * (-t3070 * t3226 + t3099 * t3377)) * t2885 * t3248) * t3071 + (t2935 * t3088 + (t2949 * t3058 - t2953 * t3062) * t3234) * t3244;
t3165 = t2885 * t3214;
t3141 = 0.2e1 * t3165;
t2815 = t3071 * t3100 * t3141 + t2827 * t3094;
t3284 = t2885 ^ 2 * t2975 * t3119;
t2820 = -t3071 * t3093 * t3284 + t2827 * t3099;
t3149 = -(t3284 + t3370) * t3071 * t3094 + t3385;
t2780 = (-t3093 * t2815 + t3099 * t3149) * t3076 + t3078 * t2820;
t3333 = t3408 * t2780;
t2823 = t2827 * t3093 + t3070 * t3284;
t2783 = (-t3099 * t2815 - t3093 * t3149) * t3076 - t3078 * t2823;
t3332 = t3408 * t2783;
t3326 = t3072 * t3095;
t3240 = t3409 * t3326;
t2866 = -pkin(5) * t2886 * t3240 + (t3072 * t3116 + t3101 * t3117) * t3366;
t2936 = -(-t3075 * t3096 + t3077 * t3312) * t3398 - pkin(5) * (t3075 * t3102 + t3077 * t3315);
t3243 = t3409 * t3338;
t3213 = t2898 * t3243;
t2828 = (-t3078 * t2866 * t3213 - (-t3095 * t3033 * t3162 + t3078 * (-t3072 * t3225 + t3101 * t3375)) * t2886 * t3243) * t3073 + (t2936 * t3088 + (t2950 * t3059 - t2954 * t3063) * t3233) * t3239;
t3163 = t2886 * t3213;
t3140 = 0.2e1 * t3163;
t2816 = t3073 * t3102 * t3140 + t2828 * t3096;
t3283 = t2886 ^ 2 * t2977 * t3119;
t2821 = -t3073 * t3095 * t3283 + t2828 * t3101;
t3148 = -(t3283 + t3369) * t3073 * t3096 + t3382;
t2781 = (-t3095 * t2816 + t3101 * t3148) * t3076 + t3078 * t2821;
t3331 = t3409 * t2781;
t2824 = t2828 * t3095 + t3072 * t3283;
t2784 = (-t3101 * t2816 - t3095 * t3148) * t3076 - t3078 * t2824;
t3330 = t3409 * t2784;
t3292 = t3236 * t3393;
t3291 = t3235 * t3390;
t3290 = t3234 * t3387;
t3289 = t3233 * t3384;
t3179 = g(1) * t3061 - g(2) * t3057;
t2998 = t3075 * t3179 + t3055;
t2852 = (t2901 * t2948 + t2902 * t2952 + t2924 * t3088) * t3407 + (t2864 * t3277 - t2869 * t3378) * t3068;
t3126 = t2852 * t3076 + t3027 * t3057 - t3028 * t3061 + t3040;
t2836 = t3092 * t2998 + t3098 * t3126;
t3288 = t2836 * t3352;
t3178 = g(1) * t3062 - g(2) * t3058;
t2999 = t3075 * t3178 + t3055;
t2853 = (t2903 * t2949 + t2904 * t2953 + t2925 * t3088) * t3408 + (t2865 * t3276 - t2870 * t3376) * t3070;
t3125 = t2853 * t3076 + t3027 * t3058 - t3028 * t3062 + t3040;
t2837 = t3094 * t2999 + t3100 * t3125;
t3287 = t2837 * t3351;
t3147 = t3019 * t3060 + t3023 * t3056;
t3274 = t3147 * t3353;
t3144 = t3022 * t3063 + t3026 * t3059;
t3273 = t3144 * t3350;
t3146 = t3020 * t3061 + t3024 * t3057;
t3272 = t3146 * t3352;
t3145 = t3021 * t3062 + t3025 * t3058;
t3271 = t3145 * t3351;
t3270 = t2931 * t3349;
t3269 = t3236 * t3147 * t3406;
t3268 = t2934 * t3346;
t3267 = t2935 * t3343;
t3266 = t2936 * t3340;
t3265 = t3235 * t3146 * t3407;
t3264 = t3234 * t3145 * t3408;
t3263 = t3233 * t3144 * t3409;
t2963 = -t3013 * t3079 + t3075 * t3324;
t3238 = t3056 * t2963 * t2831;
t3177 = g(1) * t3063 - g(2) * t3059;
t3000 = t3075 * t3177 + t3055;
t2854 = (t2905 * t2950 + t2906 * t2954 + t2926 * t3088) * t3409 + (t2866 * t3275 - t2871 * t3374) * t3072;
t3124 = t2854 * t3076 + t3027 * t3059 - t3028 * t3063 + t3040;
t2838 = t3096 * t3000 + t3102 * t3124;
t2964 = t3014 * t3095 + t3075 * t3321;
t3237 = t3059 * t2964 * t2838;
t3224 = t2809 * t3255;
t3223 = t2810 * t3250;
t3222 = t2811 * t3245;
t3221 = t2812 * t3240;
t3220 = t2831 * t3255;
t3219 = t2836 * t3261;
t3218 = t2837 * t3260;
t3217 = t2838 * t3259;
t3212 = t3066 * t3274;
t3211 = t3072 * t3273;
t3210 = t3068 * t3272;
t3209 = t3070 * t3271;
t3208 = t3236 * t3257;
t3207 = t3236 * t3256;
t3206 = t3066 * t3269;
t3205 = t3235 * t3252;
t3204 = t3235 * t3251;
t3203 = t3068 * t3265;
t3202 = t3234 * t3247;
t3201 = t3234 * t3246;
t3200 = t3070 * t3264;
t3199 = t3233 * t3242;
t3198 = t3233 * t3241;
t3197 = t3072 * t3263;
t3195 = t2836 * t3272;
t3192 = t2837 * t3271;
t3189 = t2838 * t3273;
t3187 = t2963 * t3257;
t3186 = t2964 * t3242;
t3185 = t3057 * t3288;
t3184 = t3058 * t3287;
t3183 = t3061 * t3288;
t3182 = t3062 * t3287;
t3181 = t3063 * t2838 * t3350;
t3176 = t3236 * t3224;
t3175 = t3235 * t3223;
t3174 = t3234 * t3222;
t3173 = t3233 * t3221;
t3172 = t2955 * t3220;
t3171 = t3068 * t3195;
t3170 = t3070 * t3192;
t3157 = t3072 * t3189;
t3156 = t3068 * t3185;
t3155 = t3070 * t3184;
t3154 = t3068 * t3183;
t3153 = t3070 * t3182;
t3152 = t3072 * t3181;
t3018 = -t3064 * t3087 - t3065 * t3074;
t3017 = -t3064 * t3074 + t3065 * t3087;
t2996 = t3077 * t3177 - t3397;
t2995 = t3077 * t3178 - t3397;
t2994 = t3077 * t3179 - t3397;
t2993 = t3077 * t3180 - t3397;
t2970 = t3102 * t3000;
t2969 = t3100 * t2999;
t2968 = t3098 * t2998;
t2965 = t3082 * t2997;
t2890 = (t2905 * t3026 - t2906 * t3022) * t3409;
t2889 = (t2903 * t3025 - t2904 * t3021) * t3408;
t2888 = (t2901 * t3024 - t2902 * t3020) * t3407;
t2887 = (t2899 * t3023 - t2900 * t3019) * t3406;
t2878 = (t3073 - 0.2e1) * t3369;
t2877 = (t3071 - 0.2e1) * t3370;
t2876 = (t3069 - 0.2e1) * t3371;
t2875 = (t3067 - 0.2e1) * t3373;
t2851 = -t3059 * g(1) - t3063 * g(2) + t2854;
t2850 = -t3058 * g(1) - t3062 * g(2) + t2853;
t2849 = -t3057 * g(1) - t3061 * g(2) + t2852;
t2847 = -t3056 * g(1) - t3060 * g(2) + t2848;
t2846 = -t2851 * t3078 - t2996 * t3076;
t2845 = -t2850 * t3078 - t2995 * t3076;
t2844 = -t2849 * t3078 - t2994 * t3076;
t2843 = -t2847 * t3078 - t2993 * t3076;
t2841 = t2970 + (-t2851 * t3076 + t2996 * t3078) * t3096;
t2840 = t2969 + (-t2850 * t3076 + t2995 * t3078) * t3094;
t2839 = t2968 + (-t2849 * t3076 + t2994 * t3078) * t3092;
t2835 = -t3096 * t3124 + t2970;
t2834 = -t3094 * t3125 + t2969;
t2833 = -t3092 * t3126 + t2968;
t2832 = t2965 + (-t2847 * t3076 + t2993 * t3078) * t3080;
t2829 = -t3080 * t3127 + t2965;
t2808 = t2811 * t3094 + t3100 * t3279;
t2807 = t3094 * t3279 - t3385;
t2806 = t2812 * t3096 + t3102 * t3278;
t2805 = t2810 * t3092 + t3098 * t3280;
t2804 = t3096 * t3278 - t3382;
t2803 = t3092 * t3280 - t3388;
t2802 = t2809 * t3080 + t3082 * t3282;
t2801 = t3080 * t3282 - t3391;
t2800 = (t3072 * t3140 + t3383) * t3095;
t2799 = (t3070 * t3141 + t3386) * t3093;
t2798 = (t3068 * t3142 + t3389) * t3091;
t2797 = (t3066 * t3143 + t3392) * t3079;
t2796 = t2841 * t3101 + t2846 * t3095;
t2795 = t2841 * t3095 - t2846 * t3101;
t2794 = t2840 * t3099 + t2845 * t3093;
t2793 = t2840 * t3093 - t2845 * t3099;
t2792 = t2839 * t3097 + t2844 * t3091;
t2791 = t2839 * t3091 - t2844 * t3097;
t2790 = 0.2e1 * t3101 * t3383 + (-0.2e1 * t3073 + 0.4e1) * t3163;
t2789 = 0.2e1 * t3099 * t3386 + (-0.2e1 * t3071 + 0.4e1) * t3165;
t2788 = 0.2e1 * t3097 * t3389 + (-0.2e1 * t3069 + 0.4e1) * t3167;
t2787 = 0.2e1 * t3081 * t3392 + (-0.2e1 * t3067 + 0.4e1) * t3169;
t2786 = t2832 * t3081 + t2843 * t3079;
t2785 = t2832 * t3079 - t2843 * t3081;
t1 = [(t2847 * t3364 + t2849 * t3362 + t2850 * t3360 + t2851 * t3358) * MDP(1) + (-t2809 * t3196 - t2810 * t3193 - t2811 * t3190 - t2812 * t3188) * MDP(2) + (-t2831 * t3196 - t3154 - t3153 - t3152 + (-t2801 * t3364 - t2803 * t3362 - t2804 * t3358 - t2807 * t3360) * t3076) * MDP(3) + (-t2829 * t3196 - t2833 * t3193 - t2834 * t3190 - t2835 * t3188 + (-t2802 * t3364 - t2805 * t3362 - t2806 * t3358 - t2808 * t3360) * t3076) * MDP(4) + (-t2797 * t3196 - t2798 * t3193 - t2799 * t3190 - t2800 * t3188 + (t3060 * t3135 + t3061 * t3134 + t3062 * t3133 + t3063 * t3132) * t3118) * MDP(5) + (-t2787 * t3196 - t2788 * t3193 - t2789 * t3190 - t2790 * t3188 + (-t2875 * t3207 - t2876 * t3204 - t2877 * t3201 - t2878 * t3198) * t3118) * MDP(6) + (-t2818 * t3196 - t2822 * t3193 - t2823 * t3190 - t2824 * t3188 + (-t3060 * t3176 - t3061 * t3175 - t3062 * t3174 - t3063 * t3173) * t3118) * MDP(7) + (-t2817 * t3196 - t2819 * t3193 - t2820 * t3190 - t2821 * t3188 + (-t3060 * t3292 - t3061 * t3291 - t3062 * t3290 - t3063 * t3289) * t3118) * MDP(8) + (-t2825 * t3207 - t2826 * t3204 - t2827 * t3201 - t2828 * t3198) * t3394 + (-t3060 * t2831 * t3353 - t3183 - t3182 - t3181 + t2900 * t3337 + t2902 * t3335 + t2904 * t3333 + t2906 * t3331 + (-t2785 * t3207 - t2791 * t3204 - t2793 * t3201 - t2795 * t3198) * t3118) * MDP(10) + (t2906 * t3330 + t3095 * t3152 + t2904 * t3332 + t3093 * t3153 + t2902 * t3334 + t3091 * t3154 + t2900 * t3336 + t3060 * t3172 + (-t2786 * t3207 - t2792 * t3204 - t2794 * t3201 - t2796 * t3198) * t3118) * MDP(11) + t3018 * MDP(13) - t3017 * MDP(14) + t3395 * MDP(15); (t2847 * t3365 + t2849 * t3363 + t2850 * t3361 + t2851 * t3359) * MDP(1) + (t2809 * t3187 + t2810 * t3194 + t2811 * t3191 + t2812 * t3186) * MDP(2) + (t2831 * t3187 + t3156 + t3155 + t2838 * t3186 + (-t2801 * t3365 - t2803 * t3363 - t2804 * t3359 - t2807 * t3361) * t3076) * MDP(3) + (t2829 * t3187 + t2833 * t3194 + t2834 * t3191 + t2835 * t3186 + (-t2802 * t3365 - t2805 * t3363 - t2806 * t3359 - t2808 * t3361) * t3076) * MDP(4) + (t2797 * t3187 + t2798 * t3194 + t2799 * t3191 + t2800 * t3186 + (-t3056 * t3135 - t3057 * t3134 - t3058 * t3133 - t3059 * t3132) * t3118) * MDP(5) + (t2787 * t3187 + t2788 * t3194 + t2789 * t3191 + t2790 * t3186 + (t2875 * t3208 + t2876 * t3205 + t2877 * t3202 + t2878 * t3199) * t3118) * MDP(6) + (t2818 * t3187 + t2822 * t3194 + t2823 * t3191 + t2824 * t3186 + (t3056 * t3176 + t3057 * t3175 + t3058 * t3174 + t3059 * t3173) * t3118) * MDP(7) + (t2817 * t3187 + t2819 * t3194 + t2820 * t3191 + t2821 * t3186 + (t3056 * t3292 + t3057 * t3291 + t3058 * t3290 + t3059 * t3289) * t3118) * MDP(8) + (t2825 * t3208 + t2826 * t3205 + t2827 * t3202 + t2828 * t3199) * t3394 + (t3185 + t3184 + t2901 * t3335 + t2903 * t3333 + (t2905 * t2781 + t3237) * t3409 + (t2899 * t2777 + t3238) * t3406 + (t2785 * t3208 + t2791 * t3205 + t2793 * t3202 + t2795 * t3199) * t3118) * MDP(10) + (-t3091 * t3156 - t3093 * t3155 + t2901 * t3334 + t2903 * t3332 + (t2905 * t2784 - t3237 * t3326) * t3409 + (t2899 * t2778 - t3238 * t3327) * t3406 + (t2786 * t3208 + t2792 * t3205 + t2794 * t3202 + t2796 * t3199) * t3118) * MDP(11) + t3017 * MDP(13) + t3018 * MDP(14) + t3396 * MDP(15); (t2847 * t3357 + t2849 * t3356 + t2850 * t3355 + t2851 * t3354) * MDP(1) + (t2809 * t3262 + t2810 * t3261 + t2811 * t3260 + t2812 * t3259) * MDP(2) + (t2831 * t3262 + t3219 + t3218 + t3217 + (-t2801 * t3357 - t2803 * t3356 - t2804 * t3354 - t2807 * t3355) * t3076) * MDP(3) + (t2829 * t3262 + t2833 * t3261 + t2834 * t3260 + t2835 * t3259 + (-t2802 * t3357 - t2805 * t3356 - t2806 * t3354 - t2808 * t3355) * t3076) * MDP(4) + (t2797 * t3262 + t2798 * t3261 + t2799 * t3260 + t2800 * t3259 + (-t2931 * t3161 - t2934 * t3160 - t2935 * t3159 - t2936 * t3158) * t3118) * MDP(5) + (t2787 * t3262 + t2788 * t3261 + t2789 * t3260 + t2790 * t3259 + (t2875 * t3270 + t2876 * t3268 + t2877 * t3267 + t2878 * t3266) * t3118) * MDP(6) + (t2818 * t3262 + t2822 * t3261 + t2823 * t3260 + t2824 * t3259 + (t2931 * t3224 + t2934 * t3223 + t2935 * t3222 + t2936 * t3221) * t3118) * MDP(7) + (t2817 * t3262 + t2819 * t3261 + t2820 * t3260 + t2821 * t3259 + (t2931 * t3393 + t2934 * t3390 + t2935 * t3387 + t2936 * t3384) * t3118) * MDP(8) + (t2825 * t3270 + t2826 * t3268 + t2827 * t3267 + t2828 * t3266) * t3394 + (t2920 * t3337 + t2924 * t3335 + t2925 * t3333 + t2926 * t3331 + t2956 * t3406 * t2831 + t2960 * t3407 * t2836 + t2961 * t3408 * t2837 + t2962 * t3409 * t2838 + (t2785 * t3270 + t2791 * t3268 + t2793 * t3267 + t2795 * t3266) * t3118) * MDP(10) + (t2926 * t3330 - t3095 * t3217 + t2925 * t3332 - t3093 * t3218 + t2924 * t3334 - t3091 * t3219 + t2920 * t3336 - t2956 * t3220 + (t2786 * t3270 + t2792 * t3268 + t2794 * t3267 + t2796 * t3266) * t3118) * MDP(11) + (t3088 - g(3)) * MDP(15); (t2847 * t2887 + t2849 * t2888 + t2850 * t2889 + t2851 * t2890) * MDP(1) + (t2809 * t3212 + t2810 * t3210 + t2811 * t3209 + t2812 * t3211) * MDP(2) + (t2831 * t3212 + t3171 + t3170 + t3157 + (-t2801 * t2887 - t2803 * t2888 - t2804 * t2890 - t2807 * t2889) * t3076) * MDP(3) + (t2829 * t3212 + t2833 * t3210 + t2834 * t3209 + t2835 * t3211 + (-t2802 * t2887 - t2805 * t2888 - t2806 * t2890 - t2808 * t2889) * t3076) * MDP(4) + (t2797 * t3212 + t2798 * t3210 + t2799 * t3209 + t2800 * t3211 + (-t3144 * t3132 - t3145 * t3133 - t3146 * t3134 - t3147 * t3135) * t3118) * MDP(5) + (t2787 * t3212 + t2788 * t3210 + t2789 * t3209 + t2790 * t3211 + (t2875 * t3206 + t2876 * t3203 + t2877 * t3200 + t2878 * t3197) * t3118) * MDP(6) + (t2818 * t3212 + t2822 * t3210 + t2823 * t3209 + t2824 * t3211 + (t3197 * t3383 + t3200 * t3386 + t3203 * t3389 + t3206 * t3392) * t3118) * MDP(7) + (t2817 * t3212 + t2819 * t3210 + t2820 * t3209 + t2821 * t3211 + (t2809 * t3269 + t2810 * t3265 + t2811 * t3264 + t2812 * t3263) * t3118) * MDP(8) + (t2825 * t3206 + t2826 * t3203 + t2827 * t3200 + t2828 * t3197) * t3394 + (t2831 * t3274 + t3195 + t3192 + t3189 + t2887 * t2777 + t2888 * t2779 + t2889 * t2780 + t2890 * t2781 + (t2785 * t3206 + t2791 * t3203 + t2793 * t3200 + t2795 * t3197) * t3118) * MDP(10) + (t2890 * t2784 - t3095 * t3157 + t2889 * t2783 - t3093 * t3170 + t2888 * t2782 - t3091 * t3171 + t2887 * t2778 - t3147 * t3172 + (t2786 * t3206 + t2792 * t3203 + t2794 * t3200 + t2796 * t3197) * t3118) * MDP(11) + t3087 * MDP(12) + (-t3064 * t3395 + t3065 * t3396) * MDP(13) + (-t3064 * t3396 - t3065 * t3395) * MDP(14);];
tauX  = t1;