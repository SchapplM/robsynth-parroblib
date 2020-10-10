% Calculate minimal parameter regressor of inverse dynamics forces for
% P4PRRRR8V2G3A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P4PRRRR8V2G3A0_convert_par2_MPV_fixb.m

% Output:
% tauX [4x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR8V2G3A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_mdp: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_mdp: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_mdp: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:28:32
% EndTime: 2020-08-07 11:29:03
% DurationCPUTime: 31.73s
% Computational Cost: add. (185863->871), mult. (377834->1542), div. (10692->22), fcn. (300388->30), ass. (0->569)
t3222 = sin(qJ(2,4));
t3223 = cos(qJ(3,4));
t3416 = t3222 * t3223;
t3224 = cos(qJ(2,4));
t3249 = pkin(7) + pkin(6);
t3185 = t3224 * t3249;
t3555 = pkin(2) * t3222 - t3185;
t3115 = pkin(3) * t3416 + t3555;
t3525 = pkin(3) * t3223;
t3183 = pkin(2) + t3525;
t3220 = cos(pkin(4));
t3221 = sin(qJ(3,4));
t3431 = t3220 * t3221;
t3149 = t3183 * t3431;
t3218 = sin(pkin(4));
t3444 = t3218 * t3223;
t3048 = t3115 * t3444 + t3149;
t3045 = 0.1e1 / t3048;
t3234 = sin(qJ(2,3));
t3239 = cos(qJ(3,3));
t3413 = t3234 * t3239;
t3240 = cos(qJ(2,3));
t3192 = t3240 * t3249;
t3554 = pkin(2) * t3234 - t3192;
t3116 = pkin(3) * t3413 + t3554;
t3524 = pkin(3) * t3239;
t3186 = pkin(2) + t3524;
t3233 = sin(qJ(3,3));
t3428 = t3220 * t3233;
t3150 = t3186 * t3428;
t3437 = t3218 * t3239;
t3070 = t3116 * t3437 + t3150;
t3061 = 0.1e1 / t3070;
t3236 = sin(qJ(2,2));
t3241 = cos(qJ(3,2));
t3411 = t3236 * t3241;
t3242 = cos(qJ(2,2));
t3193 = t3242 * t3249;
t3553 = pkin(2) * t3236 - t3193;
t3117 = pkin(3) * t3411 + t3553;
t3523 = pkin(3) * t3241;
t3187 = pkin(2) + t3523;
t3235 = sin(qJ(3,2));
t3426 = t3220 * t3235;
t3151 = t3187 * t3426;
t3436 = t3218 * t3241;
t3071 = t3117 * t3436 + t3151;
t3063 = 0.1e1 / t3071;
t3238 = sin(qJ(2,1));
t3243 = cos(qJ(3,1));
t3409 = t3238 * t3243;
t3244 = cos(qJ(2,1));
t3194 = t3244 * t3249;
t3552 = pkin(2) * t3238 - t3194;
t3118 = pkin(3) * t3409 + t3552;
t3522 = pkin(3) * t3243;
t3188 = pkin(2) + t3522;
t3237 = sin(qJ(3,1));
t3424 = t3220 * t3237;
t3152 = t3188 * t3424;
t3435 = t3218 * t3243;
t3072 = t3118 * t3435 + t3152;
t3065 = 0.1e1 / t3072;
t3106 = pkin(3) * t3424 + t3218 * t3552;
t3215 = t3243 ^ 2;
t3526 = pkin(3) * t3215;
t3399 = t3238 * t3526;
t3040 = 0.1e1 / (pkin(2) * t3424 + t3106 * t3243 + t3218 * t3399);
t3105 = pkin(3) * t3426 + t3218 * t3553;
t3214 = t3241 ^ 2;
t3527 = pkin(3) * t3214;
t3400 = t3236 * t3527;
t3039 = 0.1e1 / (pkin(2) * t3426 + t3105 * t3241 + t3218 * t3400);
t3104 = pkin(3) * t3428 + t3218 * t3554;
t3213 = t3239 ^ 2;
t3528 = pkin(3) * t3213;
t3401 = t3234 * t3528;
t3038 = 0.1e1 / (pkin(2) * t3428 + t3104 * t3239 + t3218 * t3401);
t3103 = pkin(3) * t3431 + t3218 * t3555;
t3212 = t3223 ^ 2;
t3529 = pkin(3) * t3212;
t3402 = t3222 * t3529;
t3037 = 0.1e1 / (pkin(2) * t3431 + t3103 * t3223 + t3218 * t3402);
t3250 = xP(4);
t3210 = sin(t3250);
t3211 = cos(t3250);
t3251 = koppelP(4,2);
t3255 = koppelP(4,1);
t3139 = t3210 * t3255 + t3211 * t3251;
t3225 = legFrame(4,2);
t3199 = sin(t3225);
t3203 = cos(t3225);
t3245 = xDP(4);
t3247 = xDP(2);
t3248 = xDP(1);
t3295 = t3210 * t3251 - t3211 * t3255;
t3024 = -(t3295 * t3245 - t3247) * t3199 + (t3139 * t3245 - t3248) * t3203;
t3217 = sin(pkin(8));
t3246 = xDP(3);
t3182 = t3246 * t3217;
t3219 = cos(pkin(8));
t3017 = t3024 * t3219 + t3182;
t3419 = t3220 * t3246;
t3430 = t3220 * t3222;
t2978 = ((t3024 * t3430 - t3224 * t3246) * t3219 + (t3024 * t3224 + t3222 * t3419) * t3217) * t3221 + t3017 * t3444;
t3507 = t2978 * t3045;
t3389 = t3249 * t3507;
t3323 = t3221 * t3389;
t3131 = t3183 * t3222 - t3185;
t3184 = t3222 * t3249;
t3432 = t3219 * t3246;
t2992 = t3017 * (t3183 * t3224 + t3184) * t3220 - (t3024 * t3217 - t3432) * t3131;
t3075 = t3131 * t3444 + t3149;
t3500 = t2992 / t3075;
t2963 = t3323 - t3500;
t3216 = t3245 ^ 2;
t3229 = xDDP(4);
t3231 = xDDP(2);
t3053 = -t3139 * t3216 - t3229 * t3295 + t3231;
t3232 = xDDP(1);
t3057 = -t3139 * t3229 + t3216 * t3295 + t3232;
t3230 = xDDP(3);
t3195 = pkin(2) ^ 2 + t3249 ^ 2;
t3259 = pkin(3) ^ 2;
t3381 = t3221 * t3500;
t3518 = 0.2e1 * pkin(2) * pkin(3);
t3327 = (-t3249 * t3381 + (t3212 * t3259 + t3223 * t3518 + t3195) * t3507) * t3037 * t3507;
t3260 = 0.1e1 / pkin(3);
t3380 = t3260 * t3500;
t3121 = t3217 * t3224 + t3219 * t3430;
t3154 = pkin(2) * t3224 + t3184;
t3446 = t3218 * t3221;
t3283 = pkin(3) * t3446 - t3220 * t3555;
t3013 = -t3121 * t3529 - t3154 * t3217 * t3223 + (pkin(2) * t3446 + t3283 * t3223) * t3219;
t3488 = t3013 * t3037;
t3041 = t3154 * t3219 + t3283 * t3217;
t3119 = t3217 * t3430 - t3219 * t3224;
t3445 = t3218 * t3222;
t3448 = t3217 * t3218;
t3534 = pkin(2) * t3221;
t3000 = (t3119 * t3199 + t3203 * t3445) * t3529 + (-t3041 * t3199 + t3103 * t3203) * t3223 + (-t3199 * t3448 + t3203 * t3220) * t3534;
t3495 = t3000 * t3037;
t2999 = (-t3119 * t3203 + t3199 * t3445) * t3529 + (t3041 * t3203 + t3103 * t3199) * t3223 + (t3199 * t3220 + t3203 * t3448) * t3534;
t3496 = t2999 * t3037;
t2940 = t3223 * t3327 + (pkin(2) * t3380 - t2963 * t3223) * t3037 * t3500 + t3053 * t3495 + t3057 * t3496 + t3230 * t3488;
t3447 = t3217 * t3220;
t3147 = g(1) * t3218 + g(2) * t3447;
t3148 = g(1) * t3447 - t3218 * g(2);
t3433 = t3219 * t3220;
t3165 = g(3) * t3433;
t3269 = t2940 * t3218 - t3147 * t3199 + t3148 * t3203 + t3165;
t2917 = t3269 * t3224;
t3315 = g(1) * t3203 - g(2) * t3199;
t3521 = t3217 * g(3);
t3099 = t3315 * t3219 - t3521;
t2915 = t3099 * t3222 + t2917;
t3252 = koppelP(3,2);
t3256 = koppelP(3,1);
t3140 = t3210 * t3256 + t3211 * t3252;
t3226 = legFrame(3,2);
t3200 = sin(t3226);
t3204 = cos(t3226);
t3294 = t3210 * t3252 - t3211 * t3256;
t3025 = -(t3294 * t3245 - t3247) * t3200 + (t3140 * t3245 - t3248) * t3204;
t3018 = t3025 * t3219 + t3182;
t3427 = t3220 * t3234;
t2982 = ((t3025 * t3427 - t3240 * t3246) * t3219 + (t3025 * t3240 + t3234 * t3419) * t3217) * t3233 + t3018 * t3437;
t3503 = t2982 * t3061;
t3387 = t3249 * t3503;
t3322 = t3233 * t3387;
t3132 = t3186 * t3234 - t3192;
t3189 = t3234 * t3249;
t2996 = (t3186 * t3240 + t3189) * t3018 * t3220 - (t3025 * t3217 - t3432) * t3132;
t3082 = t3132 * t3437 + t3150;
t3499 = t2996 / t3082;
t2965 = t3322 - t3499;
t3054 = -t3140 * t3216 - t3229 * t3294 + t3231;
t3058 = -t3140 * t3229 + t3216 * t3294 + t3232;
t3376 = t3233 * t3499;
t3326 = (-t3249 * t3376 + (t3213 * t3259 + t3239 * t3518 + t3195) * t3503) * t3038 * t3503;
t3375 = t3260 * t3499;
t3128 = t3217 * t3240 + t3219 * t3427;
t3158 = pkin(2) * t3240 + t3189;
t3443 = t3218 * t3233;
t3282 = pkin(3) * t3443 - t3220 * t3554;
t3014 = -t3128 * t3528 - t3158 * t3217 * t3239 + (pkin(2) * t3443 + t3282 * t3239) * t3219;
t3487 = t3014 * t3038;
t3042 = t3158 * t3219 + t3282 * t3217;
t3122 = t3217 * t3427 - t3219 * t3240;
t3442 = t3218 * t3234;
t3533 = pkin(2) * t3233;
t3002 = -(t3122 * t3204 - t3200 * t3442) * t3528 + (t3042 * t3204 + t3104 * t3200) * t3239 + (t3200 * t3220 + t3204 * t3448) * t3533;
t3493 = t3002 * t3038;
t3001 = (t3122 * t3200 + t3204 * t3442) * t3528 + (-t3042 * t3200 + t3104 * t3204) * t3239 + (-t3200 * t3448 + t3204 * t3220) * t3533;
t3494 = t3001 * t3038;
t2944 = t3239 * t3326 + (pkin(2) * t3375 - t2965 * t3239) * t3038 * t3499 + t3054 * t3494 + t3058 * t3493 + t3230 * t3487;
t3268 = t2944 * t3218 - t3147 * t3200 + t3148 * t3204 + t3165;
t2927 = t3268 * t3240;
t3314 = g(1) * t3204 - g(2) * t3200;
t3100 = t3314 * t3219 - t3521;
t2920 = t3100 * t3234 + t2927;
t3253 = koppelP(2,2);
t3257 = koppelP(2,1);
t3141 = t3210 * t3257 + t3211 * t3253;
t3227 = legFrame(2,2);
t3201 = sin(t3227);
t3205 = cos(t3227);
t3293 = t3210 * t3253 - t3211 * t3257;
t3026 = -(t3293 * t3245 - t3247) * t3201 + (t3141 * t3245 - t3248) * t3205;
t3019 = t3026 * t3219 + t3182;
t3425 = t3220 * t3236;
t2983 = ((t3026 * t3425 - t3242 * t3246) * t3219 + (t3026 * t3242 + t3236 * t3419) * t3217) * t3235 + t3019 * t3436;
t3502 = t2983 * t3063;
t3385 = t3249 * t3502;
t3321 = t3235 * t3385;
t3133 = t3187 * t3236 - t3193;
t3190 = t3236 * t3249;
t2997 = (t3187 * t3242 + t3190) * t3019 * t3220 - (t3026 * t3217 - t3432) * t3133;
t3083 = t3133 * t3436 + t3151;
t3498 = t2997 / t3083;
t2966 = t3321 - t3498;
t3055 = -t3141 * t3216 - t3229 * t3293 + t3231;
t3059 = -t3141 * t3229 + t3216 * t3293 + t3232;
t3374 = t3235 * t3498;
t3325 = (-t3249 * t3374 + (t3214 * t3259 + t3241 * t3518 + t3195) * t3502) * t3039 * t3502;
t3373 = t3260 * t3498;
t3129 = t3217 * t3242 + t3219 * t3425;
t3159 = pkin(2) * t3242 + t3190;
t3441 = t3218 * t3235;
t3281 = pkin(3) * t3441 - t3220 * t3553;
t3015 = -t3129 * t3527 - t3159 * t3217 * t3241 + (pkin(2) * t3441 + t3281 * t3241) * t3219;
t3486 = t3015 * t3039;
t3043 = t3159 * t3219 + t3281 * t3217;
t3123 = t3217 * t3425 - t3219 * t3242;
t3440 = t3218 * t3236;
t3532 = pkin(2) * t3235;
t3004 = -(t3123 * t3205 - t3201 * t3440) * t3527 + (t3043 * t3205 + t3105 * t3201) * t3241 + (t3201 * t3220 + t3205 * t3448) * t3532;
t3491 = t3004 * t3039;
t3003 = (t3123 * t3201 + t3205 * t3440) * t3527 + (-t3043 * t3201 + t3105 * t3205) * t3241 + (-t3201 * t3448 + t3205 * t3220) * t3532;
t3492 = t3003 * t3039;
t2945 = t3241 * t3325 + (pkin(2) * t3373 - t2966 * t3241) * t3039 * t3498 + t3055 * t3492 + t3059 * t3491 + t3230 * t3486;
t3267 = t2945 * t3218 - t3147 * t3201 + t3148 * t3205 + t3165;
t2928 = t3267 * t3242;
t3313 = g(1) * t3205 - g(2) * t3201;
t3101 = t3313 * t3219 - t3521;
t2921 = t3101 * t3236 + t2928;
t3254 = koppelP(1,2);
t3258 = koppelP(1,1);
t3142 = t3210 * t3258 + t3211 * t3254;
t3228 = legFrame(1,2);
t3202 = sin(t3228);
t3206 = cos(t3228);
t3292 = t3210 * t3254 - t3211 * t3258;
t3027 = -(t3292 * t3245 - t3247) * t3202 + (t3142 * t3245 - t3248) * t3206;
t3020 = t3027 * t3219 + t3182;
t3423 = t3220 * t3238;
t2984 = ((t3027 * t3423 - t3244 * t3246) * t3219 + (t3027 * t3244 + t3238 * t3419) * t3217) * t3237 + t3020 * t3435;
t3501 = t2984 * t3065;
t3383 = t3249 * t3501;
t3320 = t3237 * t3383;
t3134 = t3188 * t3238 - t3194;
t3191 = t3238 * t3249;
t2998 = (t3188 * t3244 + t3191) * t3020 * t3220 - (t3027 * t3217 - t3432) * t3134;
t3084 = t3134 * t3435 + t3152;
t3497 = t2998 / t3084;
t2967 = t3320 - t3497;
t3056 = -t3142 * t3216 - t3229 * t3292 + t3231;
t3060 = -t3142 * t3229 + t3216 * t3292 + t3232;
t3372 = t3237 * t3497;
t3324 = (-t3249 * t3372 + (t3215 * t3259 + t3243 * t3518 + t3195) * t3501) * t3040 * t3501;
t3371 = t3260 * t3497;
t3130 = t3217 * t3244 + t3219 * t3423;
t3160 = pkin(2) * t3244 + t3191;
t3439 = t3218 * t3237;
t3280 = pkin(3) * t3439 - t3220 * t3552;
t3016 = -t3130 * t3526 - t3160 * t3217 * t3243 + (pkin(2) * t3439 + t3280 * t3243) * t3219;
t3485 = t3016 * t3040;
t3044 = t3160 * t3219 + t3280 * t3217;
t3124 = t3217 * t3423 - t3219 * t3244;
t3438 = t3218 * t3238;
t3531 = pkin(2) * t3237;
t3006 = -(t3124 * t3206 - t3202 * t3438) * t3526 + (t3044 * t3206 + t3106 * t3202) * t3243 + (t3202 * t3220 + t3206 * t3448) * t3531;
t3489 = t3006 * t3040;
t3005 = (t3124 * t3202 + t3206 * t3438) * t3526 + (-t3044 * t3202 + t3106 * t3206) * t3243 + (-t3202 * t3448 + t3206 * t3220) * t3531;
t3490 = t3005 * t3040;
t2946 = t3243 * t3324 + (pkin(2) * t3371 - t2967 * t3243) * t3040 * t3497 + t3056 * t3490 + t3060 * t3489 + t3230 * t3485;
t3266 = t2946 * t3218 - t3147 * t3202 + t3148 * t3206 + t3165;
t2929 = t3266 * t3244;
t3312 = g(1) * t3206 - g(2) * t3202;
t3102 = t3312 * t3219 - t3521;
t2922 = t3102 * t3238 + t2929;
t3429 = t3220 * t3224;
t3021 = (t3217 * t3222 - t3219 * t3429) * t3525 + t3555 * t3217 - t3154 * t3433;
t3484 = t3021 * t3037;
t3422 = t3220 * t3240;
t3028 = (t3217 * t3234 - t3219 * t3422) * t3524 + t3554 * t3217 - t3158 * t3433;
t3414 = t3233 * t3239;
t3506 = t2982 ^ 2 / t3070 ^ 2;
t3305 = t3038 * t3414 * t3506;
t3278 = t3028 * t3305;
t3421 = t3220 * t3242;
t3030 = (t3217 * t3236 - t3219 * t3421) * t3523 + t3553 * t3217 - t3159 * t3433;
t3412 = t3235 * t3241;
t3505 = t2983 ^ 2 / t3071 ^ 2;
t3304 = t3039 * t3412 * t3505;
t3277 = t3030 * t3304;
t3420 = t3220 * t3244;
t3032 = (t3217 * t3238 - t3219 * t3420) * t3522 + t3552 * t3217 - t3160 * t3433;
t3410 = t3237 * t3243;
t3504 = t2984 ^ 2 / t3072 ^ 2;
t3303 = t3040 * t3410 * t3504;
t3276 = t3032 * t3303;
t3571 = t3218 * (t3243 * t3552 + t3399);
t3570 = t3218 * (t3241 * t3553 + t3400);
t3569 = t3218 * (t3239 * t3554 + t3401);
t3568 = t3218 * (t3223 * t3555 + t3402);
t3559 = t3183 * t3220;
t3558 = t3186 * t3220;
t3557 = t3187 * t3220;
t3556 = t3188 * t3220;
t3551 = -t3202 * t3056 + t3206 * t3060;
t3550 = -t3201 * t3055 + t3205 * t3059;
t3549 = -t3200 * t3054 + t3204 * t3058;
t3548 = -t3199 * t3053 + t3203 * t3057;
t3085 = t3119 * t3221 + t3217 * t3444;
t3086 = t3121 * t3221 + t3219 * t3444;
t3390 = t3224 * t3507;
t2895 = (-((t3218 * t3390 + t3220 * t3380) * t3529 + ((-t3381 + t3389) * t3222 + pkin(2) * t3390) * t3444 + t2963 * t3220) * t3507 + t3085 * t3230 - (t3218 * t3224 * t3380 + (t3220 * t3212 - t3416 * t3446 - t3220) * t3507) * t3500 - t3548 * t3086) * t3037;
t3543 = 0.2e1 * t2895;
t3088 = t3122 * t3233 + t3217 * t3437;
t3090 = t3128 * t3233 + t3219 * t3437;
t3388 = t3240 * t3503;
t2896 = (-((t3218 * t3388 + t3220 * t3375) * t3528 + ((-t3376 + t3387) * t3234 + pkin(2) * t3388) * t3437 + t2965 * t3220) * t3503 + t3088 * t3230 - (t3218 * t3240 * t3375 + (t3213 * t3220 - t3413 * t3443 - t3220) * t3503) * t3499 - t3549 * t3090) * t3038;
t3542 = 0.2e1 * t2896;
t3089 = t3123 * t3235 + t3217 * t3436;
t3091 = t3129 * t3235 + t3219 * t3436;
t3386 = t3242 * t3502;
t2897 = (-((t3218 * t3386 + t3220 * t3373) * t3527 + ((-t3374 + t3385) * t3236 + pkin(2) * t3386) * t3436 + t2966 * t3220) * t3502 + t3089 * t3230 - (t3218 * t3242 * t3373 + (t3214 * t3220 - t3411 * t3441 - t3220) * t3502) * t3498 - t3550 * t3091) * t3039;
t3541 = 0.2e1 * t2897;
t3087 = t3124 * t3237 + t3217 * t3435;
t3092 = t3130 * t3237 + t3219 * t3435;
t3384 = t3244 * t3501;
t2898 = (-((t3218 * t3384 + t3220 * t3371) * t3526 + ((-t3372 + t3383) * t3238 + pkin(2) * t3384) * t3435 + t2967 * t3220) * t3501 + t3087 * t3230 - (t3218 * t3244 * t3371 + (t3215 * t3220 - t3409 * t3439 - t3220) * t3501) * t3497 - t3551 * t3092) * t3040;
t3540 = 0.2e1 * t2898;
t3319 = t3045 * t3380;
t3415 = t3230 * t3260;
t3418 = t3220 * t3260;
t3434 = t3218 * t3260;
t3023 = (t3217 * t3429 + t3219 * t3222) * t3525 + t3154 * t3447 + t3219 * t3555;
t3483 = t3023 * t3037;
t3530 = pkin(2) * t3260;
t2926 = t3415 * t3483 - t3327 * t3418 - (-t3220 * t3323 + (-t3115 * t3221 * t3434 + t3220 * (t3223 * t3530 + t3212)) * t3500) * t3319 + t3548 * t3260 * t3484;
t3306 = t2978 * t3319;
t3535 = pkin(6) / 0.2e1;
t3539 = -0.2e1 * pkin(2) * t3306 - 0.2e1 * t2926 * t3535;
t3318 = t3061 * t3375;
t3034 = (t3217 * t3422 + t3219 * t3234) * t3524 + t3158 * t3447 + t3219 * t3554;
t3479 = t3034 * t3038;
t3482 = t3028 * t3038;
t2937 = t3415 * t3479 - t3326 * t3418 - (-t3220 * t3322 + (-t3116 * t3233 * t3434 + t3220 * (t3239 * t3530 + t3213)) * t3499) * t3318 + t3549 * t3260 * t3482;
t3302 = t2982 * t3318;
t3538 = -0.2e1 * pkin(2) * t3302 - 0.2e1 * t2937 * t3535;
t3317 = t3063 * t3373;
t3035 = (t3217 * t3421 + t3219 * t3236) * t3523 + t3159 * t3447 + t3219 * t3553;
t3478 = t3035 * t3039;
t3481 = t3030 * t3039;
t2938 = t3415 * t3478 - t3325 * t3418 - (-t3220 * t3321 + (-t3117 * t3235 * t3434 + t3220 * (t3241 * t3530 + t3214)) * t3498) * t3317 + t3550 * t3260 * t3481;
t3301 = t2983 * t3317;
t3537 = -0.2e1 * pkin(2) * t3301 - 0.2e1 * t2938 * t3535;
t3316 = t3065 * t3371;
t3036 = (t3217 * t3420 + t3219 * t3238) * t3522 + t3160 * t3447 + t3219 * t3552;
t3477 = t3036 * t3040;
t3480 = t3032 * t3040;
t2939 = t3415 * t3477 - t3324 * t3418 - (-t3220 * t3320 + (-t3118 * t3237 * t3434 + t3220 * (t3243 * t3530 + t3215)) * t3497) * t3316 + t3551 * t3260 * t3480;
t3300 = t2984 * t3316;
t3536 = -0.2e1 * pkin(2) * t3300 - 0.2e1 * t2939 * t3535;
t3520 = t3231 - g(2);
t3519 = t3232 - g(1);
t3517 = MDP(9) * t3260;
t3516 = t2895 * t3221;
t3515 = t2895 * t3224;
t3514 = t2896 * t3233;
t3513 = t2896 * t3240;
t3512 = t2897 * t3235;
t3511 = t2897 * t3242;
t3510 = t2898 * t3237;
t3509 = t2898 * t3244;
t3508 = t2978 ^ 2 / t3048 ^ 2;
t3049 = t3139 * t3203 - t3295 * t3199;
t3476 = t3037 * t3049;
t3475 = t3037 * t3085;
t3050 = t3140 * t3204 - t3294 * t3200;
t3474 = t3038 * t3050;
t3473 = t3038 * t3088;
t3051 = t3141 * t3205 - t3293 * t3201;
t3472 = t3039 * t3051;
t3471 = t3039 * t3089;
t3052 = t3142 * t3206 - t3292 * t3202;
t3470 = t3040 * t3052;
t3469 = t3040 * t3087;
t3468 = t3086 * t3199;
t3467 = t3086 * t3203;
t3466 = t3090 * t3200;
t3465 = t3090 * t3204;
t3464 = t3091 * t3201;
t3463 = t3091 * t3205;
t3462 = t3092 * t3202;
t3461 = t3092 * t3206;
t3460 = t3099 * t3224;
t3459 = t3100 * t3240;
t3458 = t3101 * t3242;
t3457 = t3102 * t3244;
t3417 = t3221 * t3223;
t3398 = t3037 * t3516;
t3397 = t2895 * t3037 * t3223;
t3396 = t3038 * t3514;
t3395 = t2896 * t3038 * t3239;
t3394 = t3039 * t3512;
t3393 = t2897 * t3039 * t3241;
t3392 = t3040 * t3510;
t3391 = t2898 * t3040 * t3243;
t3261 = 0.1e1 / pkin(3) ^ 2;
t3382 = t2992 ^ 2 / t3075 ^ 2 * t3261;
t3379 = t2996 ^ 2 / t3082 ^ 2 * t3261;
t3378 = t2997 ^ 2 / t3083 ^ 2 * t3261;
t3377 = t2998 ^ 2 / t3084 ^ 2 * t3261;
t3369 = t3021 * t3476;
t3368 = t3199 * t3484;
t3367 = t3203 * t3484;
t3366 = t3200 * t3482;
t3365 = t3204 * t3482;
t3363 = t3028 * t3474;
t3362 = t3201 * t3481;
t3361 = t3205 * t3481;
t3359 = t3030 * t3472;
t3358 = t3202 * t3480;
t3357 = t3206 * t3480;
t3355 = t3032 * t3470;
t3354 = t3086 * t3476;
t3353 = t3037 * t3468;
t3352 = t3037 * t3467;
t3351 = t3090 * t3474;
t3350 = t3038 * t3466;
t3349 = t3038 * t3465;
t3348 = t3091 * t3472;
t3347 = t3039 * t3464;
t3346 = t3039 * t3463;
t3345 = t3092 * t3470;
t3344 = t3040 * t3462;
t3343 = t3040 * t3461;
t3338 = t3021 * t3398;
t3337 = t3021 * t3397;
t3336 = t3028 * t3396;
t3335 = t3028 * t3395;
t3334 = t2896 * t3363;
t3333 = t3030 * t3394;
t3332 = t3030 * t3393;
t3331 = t2897 * t3359;
t3330 = t3032 * t3392;
t3329 = t3032 * t3391;
t3328 = t2898 * t3355;
t3307 = t3037 * t3417 * t3508;
t3299 = -(t3382 + t3508) * t3222 + t3515;
t3298 = -(t3379 + t3506) * t3234 + t3513;
t3297 = -(t3378 + t3505) * t3236 + t3511;
t3296 = -(t3377 + t3504) * t3238 + t3509;
t3291 = t3183 * t3448;
t3290 = t3186 * t3448;
t3289 = t3187 * t3448;
t3288 = t3188 * t3448;
t3287 = 0.2e1 * t3306;
t3286 = 0.2e1 * t3302;
t3285 = 0.2e1 * t3301;
t3284 = 0.2e1 * t3300;
t3279 = t3021 * t3307;
t3273 = pkin(2) * t3543 - pkin(6) * t3382;
t3272 = pkin(2) * t3542 - pkin(6) * t3379;
t3271 = pkin(2) * t3541 - pkin(6) * t3378;
t3270 = pkin(2) * t3540 - pkin(6) * t3377;
t3137 = pkin(2) * t3219 + t3249 * t3447;
t3138 = pkin(2) * t3447 - t3219 * t3249;
t3265 = t3119 * t3529 - (t3137 * t3224 - t3138 * t3222) * t3223;
t3264 = t3122 * t3528 - (t3137 * t3240 - t3138 * t3234) * t3239;
t3263 = t3123 * t3527 - (t3137 * t3242 - t3138 * t3236) * t3241;
t3262 = t3124 * t3526 - (t3137 * t3244 - t3138 * t3238) * t3243;
t3197 = t3219 * g(3);
t3136 = -t3210 * t3229 - t3211 * t3216;
t3135 = -t3210 * t3216 + t3211 * t3229;
t3098 = t3312 * t3217 + t3197;
t3097 = t3313 * t3217 + t3197;
t3096 = t3314 * t3217 + t3197;
t3095 = t3315 * t3217 + t3197;
t2974 = ((t3142 * t3262 - t3292 * t3571 + (-t3142 * t3288 - t3292 * t3556) * t3237) * t3206 + (-t3142 * t3571 - t3262 * t3292 + (-t3142 * t3556 + t3288 * t3292) * t3237) * t3202) * t3065;
t2973 = ((t3141 * t3263 - t3293 * t3570 + (-t3141 * t3289 - t3293 * t3557) * t3235) * t3205 + (-t3141 * t3570 - t3263 * t3293 + (-t3141 * t3557 + t3289 * t3293) * t3235) * t3201) * t3063;
t2972 = ((t3140 * t3264 - t3294 * t3569 + (-t3140 * t3290 - t3294 * t3558) * t3233) * t3204 + (-t3140 * t3569 - t3264 * t3294 + (-t3140 * t3558 + t3290 * t3294) * t3233) * t3200) * t3061;
t2971 = ((t3265 * t3139 - t3295 * t3568 + (-t3139 * t3291 - t3295 * t3559) * t3221) * t3203 + (-t3139 * t3568 - t3295 * t3265 + (-t3139 * t3559 + t3291 * t3295) * t3221) * t3199) * t3045;
t2962 = (-0.2e1 * t3215 + 0.1e1) * t3504;
t2961 = (-0.2e1 * t3214 + 0.1e1) * t3505;
t2960 = (-0.2e1 * t3213 + 0.1e1) * t3506;
t2959 = (-0.2e1 * t3212 + 0.1e1) * t3508;
t2943 = -g(1) * t3202 - g(2) * t3206 + t2946;
t2942 = -g(1) * t3201 - g(2) * t3205 + t2945;
t2941 = -g(1) * t3200 - g(2) * t3204 + t2944;
t2936 = t2943 * t3220 - t3098 * t3218;
t2935 = t2942 * t3220 - t3097 * t3218;
t2934 = t2941 * t3220 - t3096 * t3218;
t2933 = -g(1) * t3199 - g(2) * t3203 + t2940;
t2925 = -t3238 * t3266 + t3457;
t2924 = -t3236 * t3267 + t3458;
t2923 = -t3234 * t3268 + t3459;
t2919 = t2933 * t3220 - t3095 * t3218;
t2916 = -t3222 * t3269 + t3460;
t2914 = t2939 * t3243 - t3237 * t3377;
t2913 = t2939 * t3237 + t3243 * t3377;
t2912 = t2938 * t3241 - t3235 * t3378;
t2911 = t2938 * t3235 + t3241 * t3378;
t2910 = t2937 * t3239 - t3233 * t3379;
t2909 = t2937 * t3233 + t3239 * t3379;
t2908 = t2926 * t3223 - t3221 * t3382;
t2907 = t2926 * t3221 + t3223 * t3382;
t2903 = t2939 * t3238 + t3244 * t3284;
t2902 = t2938 * t3236 + t3242 * t3285;
t2901 = t2937 * t3234 + t3240 * t3286;
t2899 = t2926 * t3222 + t3224 * t3287;
t2894 = t2897 * t3236 + t3242 * t3505;
t2893 = t3236 * t3505 - t3511;
t2892 = t2898 * t3238 + t3244 * t3504;
t2891 = t2896 * t3234 + t3240 * t3506;
t2890 = t3238 * t3504 - t3509;
t2889 = t3234 * t3506 - t3513;
t2888 = t2895 * t3222 + t3224 * t3508;
t2887 = t3222 * t3508 - t3515;
t2886 = (t3243 * t3284 + t3510) * t3237;
t2885 = (t3241 * t3285 + t3512) * t3235;
t2884 = (t3239 * t3286 + t3514) * t3233;
t2883 = (t3223 * t3287 + t3516) * t3221;
t2882 = t3410 * t3540 + (0.4e1 * t3215 - 0.2e1) * t3300;
t2881 = t3412 * t3541 + (0.4e1 * t3214 - 0.2e1) * t3301;
t2880 = t3414 * t3542 + (0.4e1 * t3213 - 0.2e1) * t3302;
t2879 = t3417 * t3543 + (0.4e1 * t3212 - 0.2e1) * t3306;
t2878 = pkin(2) * t3504 - pkin(6) * t2898 + (-t2943 * t3218 - t3098 * t3220) * t3238 + t3457;
t2877 = pkin(2) * t3505 - pkin(6) * t2897 + (-t2942 * t3218 - t3097 * t3220) * t3236 + t3458;
t2876 = pkin(2) * t3506 - pkin(6) * t2896 + (-t2941 * t3218 - t3096 * t3220) * t3234 + t3459;
t2875 = pkin(2) * t3508 - pkin(6) * t2895 + (-t2933 * t3218 - t3095 * t3220) * t3222 + t3460;
t2874 = t2878 * t3243 - t2936 * t3237;
t2873 = t2878 * t3237 + t2936 * t3243;
t2872 = t2877 * t3241 - t2935 * t3235;
t2871 = t2877 * t3235 + t2935 * t3241;
t2870 = t2876 * t3239 - t2934 * t3233;
t2869 = t2876 * t3233 + t2934 * t3239;
t2868 = t2875 * t3223 - t2919 * t3221;
t2867 = t2875 * t3221 + t2919 * t3223;
t2866 = t3243 * t3536 + (-t3270 - t2922) * t3237;
t2865 = t3241 * t3537 + (-t3271 - t2921) * t3235;
t2864 = t3239 * t3538 + (-t3272 - t2920) * t3233;
t2863 = (t2929 + t3270) * t3243 + t3237 * t3536 + t3102 * t3409;
t2862 = (t2928 + t3271) * t3241 + t3235 * t3537 + t3101 * t3411;
t2861 = (t2927 + t3272) * t3239 + t3233 * t3538 + t3100 * t3413;
t2860 = t3223 * t3539 + (-t3273 - t2915) * t3221;
t2859 = (t2917 + t3273) * t3223 + t3221 * t3539 + t3099 * t3416;
t2858 = (-t3237 * t2903 + t3296 * t3243) * t3218 + t3220 * t2914;
t2857 = (-t3243 * t2903 - t3296 * t3237) * t3218 - t3220 * t2913;
t2856 = (-t3235 * t2902 + t3297 * t3241) * t3218 + t3220 * t2912;
t2855 = (-t3241 * t2902 - t3297 * t3235) * t3218 - t3220 * t2911;
t2854 = (-t3233 * t2901 + t3298 * t3239) * t3218 + t3220 * t2910;
t2853 = (-t3239 * t2901 - t3298 * t3233) * t3218 - t3220 * t2909;
t2852 = (-t3221 * t2899 + t3299 * t3223) * t3218 + t3220 * t2908;
t2851 = (-t3223 * t2899 - t3299 * t3221) * t3218 - t3220 * t2907;
t1 = [(t2933 * t3496 + t2941 * t3493 + t2942 * t3491 + t2943 * t3489) * MDP(1) + (-t2895 * t3352 - t2896 * t3349 - t2897 * t3346 - t2898 * t3343) * MDP(2) + (-t2915 * t3352 - t2920 * t3349 - t2921 * t3346 - t2922 * t3343 + (-t2887 * t3496 - t2889 * t3493 - t2890 * t3489 - t2893 * t3491) * t3218) * MDP(3) + (-t2916 * t3352 - t2923 * t3349 - t2924 * t3346 - t2925 * t3343 + (-t2888 * t3496 - t2891 * t3493 - t2892 * t3489 - t2894 * t3491) * t3218) * MDP(4) + (-t2883 * t3352 - t2884 * t3349 - t2885 * t3346 - t2886 * t3343 + (-t3203 * t3279 - t3204 * t3278 - t3205 * t3277 - t3206 * t3276) * t3260) * MDP(5) + (-t2879 * t3352 - t2880 * t3349 - t2881 * t3346 - t2882 * t3343 + (t2959 * t3367 + t2960 * t3365 + t2961 * t3361 + t2962 * t3357) * t3260) * MDP(6) + (-t2907 * t3352 - t2909 * t3349 - t2911 * t3346 - t2913 * t3343 + (t3203 * t3338 + t3204 * t3336 + t3205 * t3333 + t3206 * t3330) * t3260) * MDP(7) + (-t2908 * t3352 - t2910 * t3349 - t2912 * t3346 - t2914 * t3343 + (t3203 * t3337 + t3204 * t3335 + t3205 * t3332 + t3206 * t3329) * t3260) * MDP(8) + (t2926 * t3367 + t2937 * t3365 + t2938 * t3361 + t2939 * t3357) * t3517 + ((t3006 * t2858 - t2863 * t3461) * t3040 + (t3004 * t2856 - t2862 * t3463) * t3039 + (t3002 * t2854 - t2861 * t3465) * t3038 + (t2999 * t2852 - t2859 * t3467) * t3037 + (t2867 * t3367 + t2869 * t3365 + t2871 * t3361 + t2873 * t3357) * t3260) * MDP(10) + ((t3006 * t2857 - t2866 * t3461) * t3040 + (t3004 * t2855 - t2865 * t3463) * t3039 + (t3002 * t2853 - t2864 * t3465) * t3038 + (t2999 * t2851 - t2860 * t3467) * t3037 + (t2868 * t3367 + t2870 * t3365 + t2872 * t3361 + t2874 * t3357) * t3260) * MDP(11) + t3136 * MDP(13) - t3135 * MDP(14) + t3519 * MDP(15); (t2933 * t3495 + t2941 * t3494 + t2942 * t3492 + t2943 * t3490) * MDP(1) + (t2895 * t3353 + t2896 * t3350 + t2897 * t3347 + t2898 * t3344) * MDP(2) + (t2915 * t3353 + t2920 * t3350 + t2921 * t3347 + t2922 * t3344 + (-t2887 * t3495 - t2889 * t3494 - t2890 * t3490 - t2893 * t3492) * t3218) * MDP(3) + (t2916 * t3353 + t2923 * t3350 + t2924 * t3347 + t2925 * t3344 + (-t2888 * t3495 - t2891 * t3494 - t2892 * t3490 - t2894 * t3492) * t3218) * MDP(4) + (t2883 * t3353 + t2884 * t3350 + t2885 * t3347 + t2886 * t3344 + (t3199 * t3279 + t3200 * t3278 + t3201 * t3277 + t3202 * t3276) * t3260) * MDP(5) + (t2879 * t3353 + t2880 * t3350 + t2881 * t3347 + t2882 * t3344 + (-t2959 * t3368 - t2960 * t3366 - t2961 * t3362 - t2962 * t3358) * t3260) * MDP(6) + (t2907 * t3353 + t2909 * t3350 + t2911 * t3347 + t2913 * t3344 + (-t3199 * t3338 - t3200 * t3336 - t3201 * t3333 - t3202 * t3330) * t3260) * MDP(7) + (t2908 * t3353 + t2910 * t3350 + t2912 * t3347 + t2914 * t3344 + (-t3199 * t3337 - t3200 * t3335 - t3201 * t3332 - t3202 * t3329) * t3260) * MDP(8) + (-t2926 * t3368 - t2937 * t3366 - t2938 * t3362 - t2939 * t3358) * t3517 + ((t3005 * t2858 + t2863 * t3462) * t3040 + (t3003 * t2856 + t2862 * t3464) * t3039 + (t3001 * t2854 + t2861 * t3466) * t3038 + (t3000 * t2852 + t2859 * t3468) * t3037 + (-t2867 * t3368 - t2869 * t3366 - t2871 * t3362 - t2873 * t3358) * t3260) * MDP(10) + ((t3005 * t2857 + t2866 * t3462) * t3040 + (t3003 * t2855 + t2865 * t3464) * t3039 + (t3001 * t2853 + t2864 * t3466) * t3038 + (t3000 * t2851 + t2860 * t3468) * t3037 + (-t2868 * t3368 - t2870 * t3366 - t2872 * t3362 - t2874 * t3358) * t3260) * MDP(11) + t3135 * MDP(13) + t3136 * MDP(14) + t3520 * MDP(15); (t2933 * t3488 + t2941 * t3487 + t2942 * t3486 + t2943 * t3485) * MDP(1) + (t2895 * t3475 + t2896 * t3473 + t2897 * t3471 + t2898 * t3469) * MDP(2) + (t2915 * t3475 + t2920 * t3473 + t2921 * t3471 + t2922 * t3469 + (-t2887 * t3488 - t2889 * t3487 - t2890 * t3485 - t2893 * t3486) * t3218) * MDP(3) + (t2916 * t3475 + t2923 * t3473 + t2924 * t3471 + t2925 * t3469 + (-t2888 * t3488 - t2891 * t3487 - t2892 * t3485 - t2894 * t3486) * t3218) * MDP(4) + (t2883 * t3475 + t2884 * t3473 + t2885 * t3471 + t2886 * t3469 + (-t3023 * t3307 - t3034 * t3305 - t3035 * t3304 - t3036 * t3303) * t3260) * MDP(5) + (t2879 * t3475 + t2880 * t3473 + t2881 * t3471 + t2882 * t3469 + (t2959 * t3483 + t2960 * t3479 + t2961 * t3478 + t2962 * t3477) * t3260) * MDP(6) + (t2907 * t3475 + t2909 * t3473 + t2911 * t3471 + t2913 * t3469 + (t3023 * t3398 + t3034 * t3396 + t3035 * t3394 + t3036 * t3392) * t3260) * MDP(7) + (t2908 * t3475 + t2910 * t3473 + t2912 * t3471 + t2914 * t3469 + (t3023 * t3397 + t3034 * t3395 + t3035 * t3393 + t3036 * t3391) * t3260) * MDP(8) + (t2926 * t3483 + t2937 * t3479 + t2938 * t3478 + t2939 * t3477) * t3517 + ((t2858 * t3016 + t2863 * t3087) * t3040 + (t2856 * t3015 + t2862 * t3089) * t3039 + (t2854 * t3014 + t2861 * t3088) * t3038 + (t2852 * t3013 + t2859 * t3085) * t3037 + (t2867 * t3483 + t2869 * t3479 + t2871 * t3478 + t2873 * t3477) * t3260) * MDP(10) + ((t2857 * t3016 + t2866 * t3087) * t3040 + (t2855 * t3015 + t2865 * t3089) * t3039 + (t2853 * t3014 + t2864 * t3088) * t3038 + (t2851 * t3013 + t2860 * t3085) * t3037 + (t2868 * t3483 + t2870 * t3479 + t2872 * t3478 + t2874 * t3477) * t3260) * MDP(11) + (t3230 - g(3)) * MDP(15); (t2933 * t2971 + t2941 * t2972 + t2942 * t2973 + t2943 * t2974) * MDP(1) + (t2895 * t3354 + t2896 * t3351 + t2897 * t3348 + t2898 * t3345) * MDP(2) + (t2915 * t3354 + t2920 * t3351 + t2921 * t3348 + t2922 * t3345 + (-t2887 * t2971 - t2889 * t2972 - t2890 * t2974 - t2893 * t2973) * t3218) * MDP(3) + (t2916 * t3354 + t2923 * t3351 + t2924 * t3348 + t2925 * t3345 + (-t2888 * t2971 - t2891 * t2972 - t2892 * t2974 - t2894 * t2973) * t3218) * MDP(4) + (t2883 * t3354 + t2884 * t3351 + t2885 * t3348 + t2886 * t3345 + (t3049 * t3279 + t3050 * t3278 + t3051 * t3277 + t3052 * t3276) * t3260) * MDP(5) + (t2879 * t3354 + t2880 * t3351 + t2881 * t3348 + t2882 * t3345 + (-t2959 * t3369 - t2960 * t3363 - t2961 * t3359 - t2962 * t3355) * t3260) * MDP(6) + (t2907 * t3354 + t2909 * t3351 + t2911 * t3348 + t2913 * t3345 + (-t3049 * t3338 - t3233 * t3334 - t3235 * t3331 - t3237 * t3328) * t3260) * MDP(7) + (t2908 * t3354 + t2910 * t3351 + t2912 * t3348 + t2914 * t3345 + (-t3049 * t3337 - t3239 * t3334 - t3241 * t3331 - t3243 * t3328) * t3260) * MDP(8) + (-t2926 * t3369 - t2937 * t3363 - t2938 * t3359 - t2939 * t3355) * t3517 + (t2859 * t3354 + t2861 * t3351 + t2862 * t3348 + t2863 * t3345 + t2971 * t2852 + t2972 * t2854 + t2973 * t2856 + t2974 * t2858 + (-t2867 * t3369 - t2869 * t3363 - t2871 * t3359 - t2873 * t3355) * t3260) * MDP(10) + (t2860 * t3354 + t2864 * t3351 + t2865 * t3348 + t2866 * t3345 + t2971 * t2851 + t2972 * t2853 + t2973 * t2855 + t2974 * t2857 + (-t2868 * t3369 - t2870 * t3363 - t2872 * t3359 - t2874 * t3355) * t3260) * MDP(11) + t3229 * MDP(12) + (-t3519 * t3210 + t3520 * t3211) * MDP(13) + (-t3520 * t3210 - t3519 * t3211) * MDP(14);];
tauX  = t1;
