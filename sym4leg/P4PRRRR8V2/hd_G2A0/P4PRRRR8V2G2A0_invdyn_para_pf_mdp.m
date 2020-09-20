% Calculate minimal parameter regressor of inverse dynamics forces for
% P4PRRRR8V2G2A0
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
%   see P4PRRRR8V2G2A0_convert_par2_MPV_fixb.m

% Output:
% tauX [4x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:22
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR8V2G2A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_mdp: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_mdp: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_mdp: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:19:47
% EndTime: 2020-08-07 11:20:17
% DurationCPUTime: 31.54s
% Computational Cost: add. (185863->871), mult. (377834->1539), div. (10692->22), fcn. (300388->30), ass. (0->570)
t3237 = sin(qJ(2,4));
t3238 = cos(qJ(3,4));
t3433 = t3237 * t3238;
t3239 = cos(qJ(2,4));
t3264 = pkin(7) + pkin(6);
t3200 = t3239 * t3264;
t3576 = pkin(2) * t3237 - t3200;
t3138 = pkin(3) * t3433 + t3576;
t3542 = pkin(3) * t3238;
t3198 = pkin(2) + t3542;
t3235 = cos(pkin(4));
t3236 = sin(qJ(3,4));
t3448 = t3235 * t3236;
t3168 = t3198 * t3448;
t3233 = sin(pkin(4));
t3465 = t3233 * t3238;
t3069 = t3138 * t3465 + t3168;
t3066 = 0.1e1 / t3069;
t3249 = sin(qJ(2,3));
t3254 = cos(qJ(3,3));
t3430 = t3249 * t3254;
t3255 = cos(qJ(2,3));
t3207 = t3255 * t3264;
t3575 = pkin(2) * t3249 - t3207;
t3139 = pkin(3) * t3430 + t3575;
t3541 = pkin(3) * t3254;
t3201 = pkin(2) + t3541;
t3248 = sin(qJ(3,3));
t3445 = t3235 * t3248;
t3169 = t3201 * t3445;
t3458 = t3233 * t3254;
t3091 = t3139 * t3458 + t3169;
t3082 = 0.1e1 / t3091;
t3251 = sin(qJ(2,2));
t3256 = cos(qJ(3,2));
t3428 = t3251 * t3256;
t3257 = cos(qJ(2,2));
t3208 = t3257 * t3264;
t3574 = pkin(2) * t3251 - t3208;
t3140 = pkin(3) * t3428 + t3574;
t3540 = pkin(3) * t3256;
t3202 = pkin(2) + t3540;
t3250 = sin(qJ(3,2));
t3443 = t3235 * t3250;
t3170 = t3202 * t3443;
t3457 = t3233 * t3256;
t3092 = t3140 * t3457 + t3170;
t3084 = 0.1e1 / t3092;
t3253 = sin(qJ(2,1));
t3258 = cos(qJ(3,1));
t3426 = t3253 * t3258;
t3259 = cos(qJ(2,1));
t3209 = t3259 * t3264;
t3573 = pkin(2) * t3253 - t3209;
t3141 = pkin(3) * t3426 + t3573;
t3539 = pkin(3) * t3258;
t3203 = pkin(2) + t3539;
t3252 = sin(qJ(3,1));
t3441 = t3235 * t3252;
t3171 = t3203 * t3441;
t3456 = t3233 * t3258;
t3093 = t3141 * t3456 + t3171;
t3086 = 0.1e1 / t3093;
t3129 = pkin(3) * t3441 + t3233 * t3573;
t3230 = t3258 ^ 2;
t3543 = pkin(3) * t3230;
t3418 = t3253 * t3543;
t3061 = 0.1e1 / (pkin(2) * t3441 + t3129 * t3258 + t3233 * t3418);
t3128 = pkin(3) * t3443 + t3233 * t3574;
t3229 = t3256 ^ 2;
t3544 = pkin(3) * t3229;
t3419 = t3251 * t3544;
t3060 = 0.1e1 / (pkin(2) * t3443 + t3128 * t3256 + t3233 * t3419);
t3127 = pkin(3) * t3445 + t3233 * t3575;
t3228 = t3254 ^ 2;
t3545 = pkin(3) * t3228;
t3420 = t3249 * t3545;
t3059 = 0.1e1 / (pkin(2) * t3445 + t3127 * t3254 + t3233 * t3420);
t3126 = pkin(3) * t3448 + t3233 * t3576;
t3227 = t3238 ^ 2;
t3546 = pkin(3) * t3227;
t3421 = t3237 * t3546;
t3058 = 0.1e1 / (pkin(2) * t3448 + t3126 * t3238 + t3233 * t3421);
t3265 = xP(4);
t3225 = sin(t3265);
t3226 = cos(t3265);
t3266 = koppelP(4,2);
t3270 = koppelP(4,1);
t3158 = t3225 * t3270 + t3226 * t3266;
t3240 = legFrame(4,2);
t3214 = sin(t3240);
t3218 = cos(t3240);
t3260 = xDP(4);
t3262 = xDP(2);
t3263 = xDP(1);
t3317 = t3225 * t3266 - t3226 * t3270;
t3048 = -(t3317 * t3260 - t3262) * t3214 + (t3158 * t3260 - t3263) * t3218;
t3232 = sin(pkin(8));
t3234 = cos(pkin(8));
t3261 = xDP(3);
t3449 = t3234 * t3261;
t3042 = t3048 * t3232 - t3449;
t3436 = t3235 * t3261;
t3447 = t3235 * t3237;
t3005 = ((t3048 * t3447 - t3239 * t3261) * t3232 - (t3048 * t3239 + t3237 * t3436) * t3234) * t3236 + t3042 * t3465;
t3524 = t3005 * t3066;
t3408 = t3264 * t3524;
t3345 = t3236 * t3408;
t3150 = t3198 * t3237 - t3200;
t3197 = t3261 * t3232;
t3199 = t3237 * t3264;
t3019 = (t3198 * t3239 + t3199) * t3042 * t3235 + t3150 * (t3048 * t3234 + t3197);
t3096 = t3150 * t3465 + t3168;
t3517 = t3019 / t3096;
t2990 = t3345 - t3517;
t3231 = t3260 ^ 2;
t3244 = xDDP(4);
t3246 = xDDP(2);
t3074 = -t3158 * t3231 - t3244 * t3317 + t3246;
t3247 = xDDP(1);
t3078 = -t3158 * t3244 + t3231 * t3317 + t3247;
t3245 = xDDP(3);
t3210 = pkin(2) ^ 2 + t3264 ^ 2;
t3274 = pkin(3) ^ 2;
t3400 = t3236 * t3517;
t3535 = 0.2e1 * pkin(2) * pkin(3);
t3349 = (-t3264 * t3400 + (t3227 * t3274 + t3238 * t3535 + t3210) * t3524) * t3058 * t3524;
t3275 = 0.1e1 / pkin(3);
t3399 = t3275 * t3517;
t3142 = t3232 * t3447 - t3234 * t3239;
t3173 = pkin(2) * t3239 + t3199;
t3467 = t3233 * t3236;
t3305 = pkin(3) * t3467 - t3235 * t3576;
t3453 = t3234 * t3238;
t3038 = -t3142 * t3546 + t3173 * t3453 + (pkin(2) * t3467 + t3305 * t3238) * t3232;
t3505 = t3038 * t3058;
t3062 = -t3173 * t3232 + t3305 * t3234;
t3143 = t3232 * t3239 + t3234 * t3447;
t3466 = t3233 * t3237;
t3468 = t3233 * t3234;
t3551 = pkin(2) * t3236;
t3027 = -(t3143 * t3214 - t3218 * t3466) * t3546 + (t3062 * t3214 + t3126 * t3218) * t3238 + (t3214 * t3468 + t3218 * t3235) * t3551;
t3512 = t3027 * t3058;
t3026 = (t3143 * t3218 + t3214 * t3466) * t3546 + (-t3062 * t3218 + t3126 * t3214) * t3238 + (t3214 * t3235 - t3218 * t3468) * t3551;
t3513 = t3026 * t3058;
t2967 = t3238 * t3349 + (pkin(2) * t3399 - t2990 * t3238) * t3058 * t3517 + t3078 * t3513 + t3074 * t3512 + t3245 * t3505;
t3454 = t3234 * t3235;
t3166 = -t3233 * g(1) + g(2) * t3454;
t3167 = g(1) * t3454 + g(2) * t3233;
t3469 = t3232 * t3235;
t3184 = g(3) * t3469;
t3284 = t2967 * t3233 + t3166 * t3214 - t3167 * t3218 + t3184;
t2944 = t3284 * t3239;
t3212 = g(3) * t3234;
t3337 = g(1) * t3218 - g(2) * t3214;
t3569 = t3337 * t3232 + t3212;
t2943 = t3237 * t3569 + t2944;
t3267 = koppelP(3,2);
t3271 = koppelP(3,1);
t3159 = t3225 * t3271 + t3226 * t3267;
t3241 = legFrame(3,2);
t3215 = sin(t3241);
t3219 = cos(t3241);
t3316 = t3225 * t3267 - t3226 * t3271;
t3049 = -(t3316 * t3260 - t3262) * t3215 + (t3159 * t3260 - t3263) * t3219;
t3043 = t3049 * t3232 - t3449;
t3444 = t3235 * t3249;
t3009 = ((t3049 * t3444 - t3255 * t3261) * t3232 - (t3049 * t3255 + t3249 * t3436) * t3234) * t3248 + t3043 * t3458;
t3520 = t3009 * t3082;
t3406 = t3264 * t3520;
t3344 = t3248 * t3406;
t3151 = t3201 * t3249 - t3207;
t3204 = t3249 * t3264;
t3023 = (t3201 * t3255 + t3204) * t3043 * t3235 + (t3049 * t3234 + t3197) * t3151;
t3103 = t3151 * t3458 + t3169;
t3516 = t3023 / t3103;
t2992 = t3344 - t3516;
t3075 = -t3159 * t3231 - t3244 * t3316 + t3246;
t3079 = -t3159 * t3244 + t3231 * t3316 + t3247;
t3395 = t3248 * t3516;
t3348 = (-t3264 * t3395 + (t3228 * t3274 + t3254 * t3535 + t3210) * t3520) * t3059 * t3520;
t3394 = t3275 * t3516;
t3144 = t3232 * t3444 - t3234 * t3255;
t3177 = pkin(2) * t3255 + t3204;
t3464 = t3233 * t3248;
t3304 = pkin(3) * t3464 - t3235 * t3575;
t3452 = t3234 * t3254;
t3039 = -t3144 * t3545 + t3177 * t3452 + (pkin(2) * t3464 + t3304 * t3254) * t3232;
t3504 = t3039 * t3059;
t3063 = -t3177 * t3232 + t3304 * t3234;
t3147 = t3232 * t3255 + t3234 * t3444;
t3463 = t3233 * t3249;
t3550 = pkin(2) * t3248;
t3031 = -(t3147 * t3215 - t3219 * t3463) * t3545 + (t3063 * t3215 + t3127 * t3219) * t3254 + (t3215 * t3468 + t3219 * t3235) * t3550;
t3508 = t3031 * t3059;
t3028 = (t3147 * t3219 + t3215 * t3463) * t3545 + (-t3063 * t3219 + t3127 * t3215) * t3254 + (t3215 * t3235 - t3219 * t3468) * t3550;
t3511 = t3028 * t3059;
t2971 = t3254 * t3348 + (pkin(2) * t3394 - t2992 * t3254) * t3059 * t3516 + t3079 * t3511 + t3075 * t3508 + t3245 * t3504;
t3283 = t2971 * t3233 + t3166 * t3215 - t3167 * t3219 + t3184;
t2954 = t3283 * t3255;
t3336 = g(1) * t3219 - g(2) * t3215;
t3570 = t3336 * t3232 + t3212;
t2950 = t3249 * t3570 + t2954;
t3268 = koppelP(2,2);
t3272 = koppelP(2,1);
t3160 = t3225 * t3272 + t3226 * t3268;
t3242 = legFrame(2,2);
t3216 = sin(t3242);
t3220 = cos(t3242);
t3315 = t3225 * t3268 - t3226 * t3272;
t3050 = -(t3315 * t3260 - t3262) * t3216 + (t3160 * t3260 - t3263) * t3220;
t3044 = t3050 * t3232 - t3449;
t3442 = t3235 * t3251;
t3010 = ((t3050 * t3442 - t3257 * t3261) * t3232 - (t3050 * t3257 + t3251 * t3436) * t3234) * t3250 + t3044 * t3457;
t3519 = t3010 * t3084;
t3404 = t3264 * t3519;
t3343 = t3250 * t3404;
t3152 = t3202 * t3251 - t3208;
t3205 = t3251 * t3264;
t3024 = (t3202 * t3257 + t3205) * t3044 * t3235 + (t3050 * t3234 + t3197) * t3152;
t3104 = t3152 * t3457 + t3170;
t3515 = t3024 / t3104;
t2993 = t3343 - t3515;
t3076 = -t3160 * t3231 - t3244 * t3315 + t3246;
t3080 = -t3160 * t3244 + t3231 * t3315 + t3247;
t3393 = t3250 * t3515;
t3347 = (-t3264 * t3393 + (t3229 * t3274 + t3256 * t3535 + t3210) * t3519) * t3060 * t3519;
t3392 = t3275 * t3515;
t3145 = t3232 * t3442 - t3234 * t3257;
t3178 = pkin(2) * t3257 + t3205;
t3462 = t3233 * t3250;
t3303 = pkin(3) * t3462 - t3235 * t3574;
t3451 = t3234 * t3256;
t3040 = -t3145 * t3544 + t3178 * t3451 + (pkin(2) * t3462 + t3303 * t3256) * t3232;
t3503 = t3040 * t3060;
t3064 = -t3178 * t3232 + t3303 * t3234;
t3148 = t3232 * t3257 + t3234 * t3442;
t3461 = t3233 * t3251;
t3549 = pkin(2) * t3250;
t3032 = -(t3148 * t3216 - t3220 * t3461) * t3544 + (t3064 * t3216 + t3128 * t3220) * t3256 + (t3216 * t3468 + t3220 * t3235) * t3549;
t3507 = t3032 * t3060;
t3029 = (t3148 * t3220 + t3216 * t3461) * t3544 + (-t3064 * t3220 + t3128 * t3216) * t3256 + (t3216 * t3235 - t3220 * t3468) * t3549;
t3510 = t3029 * t3060;
t2972 = t3256 * t3347 + (pkin(2) * t3392 - t2993 * t3256) * t3060 * t3515 + t3080 * t3510 + t3076 * t3507 + t3245 * t3503;
t3282 = t2972 * t3233 + t3166 * t3216 - t3167 * t3220 + t3184;
t2955 = t3282 * t3257;
t3335 = g(1) * t3220 - g(2) * t3216;
t3571 = t3335 * t3232 + t3212;
t2951 = t3251 * t3571 + t2955;
t3269 = koppelP(1,2);
t3273 = koppelP(1,1);
t3161 = t3225 * t3273 + t3226 * t3269;
t3243 = legFrame(1,2);
t3217 = sin(t3243);
t3221 = cos(t3243);
t3314 = t3225 * t3269 - t3226 * t3273;
t3051 = -(t3314 * t3260 - t3262) * t3217 + (t3161 * t3260 - t3263) * t3221;
t3045 = t3051 * t3232 - t3449;
t3440 = t3235 * t3253;
t3011 = ((t3051 * t3440 - t3259 * t3261) * t3232 - (t3051 * t3259 + t3253 * t3436) * t3234) * t3252 + t3045 * t3456;
t3518 = t3011 * t3086;
t3402 = t3264 * t3518;
t3342 = t3252 * t3402;
t3153 = t3203 * t3253 - t3209;
t3206 = t3253 * t3264;
t3025 = (t3203 * t3259 + t3206) * t3045 * t3235 + (t3051 * t3234 + t3197) * t3153;
t3105 = t3153 * t3456 + t3171;
t3514 = t3025 / t3105;
t2994 = t3342 - t3514;
t3077 = -t3161 * t3231 - t3244 * t3314 + t3246;
t3081 = -t3161 * t3244 + t3231 * t3314 + t3247;
t3391 = t3252 * t3514;
t3346 = (-t3264 * t3391 + (t3230 * t3274 + t3258 * t3535 + t3210) * t3518) * t3061 * t3518;
t3390 = t3275 * t3514;
t3146 = t3232 * t3440 - t3234 * t3259;
t3179 = pkin(2) * t3259 + t3206;
t3460 = t3233 * t3252;
t3302 = pkin(3) * t3460 - t3235 * t3573;
t3450 = t3234 * t3258;
t3041 = -t3146 * t3543 + t3179 * t3450 + (pkin(2) * t3460 + t3302 * t3258) * t3232;
t3502 = t3041 * t3061;
t3065 = -t3179 * t3232 + t3302 * t3234;
t3149 = t3232 * t3259 + t3234 * t3440;
t3459 = t3233 * t3253;
t3548 = pkin(2) * t3252;
t3033 = -(t3149 * t3217 - t3221 * t3459) * t3543 + (t3065 * t3217 + t3129 * t3221) * t3258 + (t3217 * t3468 + t3221 * t3235) * t3548;
t3506 = t3033 * t3061;
t3030 = (t3149 * t3221 + t3217 * t3459) * t3543 + (-t3065 * t3221 + t3129 * t3217) * t3258 + (t3217 * t3235 - t3221 * t3468) * t3548;
t3509 = t3030 * t3061;
t2973 = t3258 * t3346 + (pkin(2) * t3390 - t2994 * t3258) * t3061 * t3514 + t3081 * t3509 + t3077 * t3506 + t3245 * t3502;
t3281 = t2973 * t3233 + t3166 * t3217 - t3167 * t3221 + t3184;
t2956 = t3281 * t3259;
t3334 = g(1) * t3221 - g(2) * t3217;
t3572 = t3334 * t3232 + t3212;
t2952 = t3253 * t3572 + t2956;
t3592 = (t3258 * t3573 + t3418) * t3233;
t3591 = (t3256 * t3574 + t3419) * t3233;
t3590 = (t3254 * t3575 + t3420) * t3233;
t3589 = (t3238 * t3576 + t3421) * t3233;
t3588 = t3569 * t3239;
t3587 = t3570 * t3255;
t3586 = t3571 * t3257;
t3585 = t3572 * t3259;
t3580 = t3198 * t3235;
t3579 = t3201 * t3235;
t3578 = t3202 * t3235;
t3577 = t3203 * t3235;
t3568 = t3077 * t3217 - t3081 * t3221;
t3567 = t3076 * t3216 - t3080 * t3220;
t3566 = t3075 * t3215 - t3079 * t3219;
t3565 = t3074 * t3214 - t3078 * t3218;
t3108 = -t3143 * t3236 - t3233 * t3453;
t3293 = t3142 * t3236 + t3232 * t3465;
t3409 = t3239 * t3524;
t2922 = (-((t3233 * t3409 + t3235 * t3399) * t3546 + ((-t3400 + t3408) * t3237 + pkin(2) * t3409) * t3465 + t2990 * t3235) * t3524 + t3108 * t3245 - (t3233 * t3239 * t3399 + (t3227 * t3235 - t3433 * t3467 - t3235) * t3524) * t3517 + t3565 * t3293) * t3058;
t3560 = 0.2e1 * t2922;
t3109 = t3144 * t3248 + t3232 * t3458;
t3115 = -t3147 * t3248 - t3233 * t3452;
t3407 = t3255 * t3520;
t2923 = (-((t3233 * t3407 + t3235 * t3394) * t3545 + ((-t3395 + t3406) * t3249 + pkin(2) * t3407) * t3458 + t2992 * t3235) * t3520 + t3115 * t3245 - (t3233 * t3255 * t3394 + (t3228 * t3235 - t3430 * t3464 - t3235) * t3520) * t3516 + t3566 * t3109) * t3059;
t3559 = 0.2e1 * t2923;
t3111 = t3145 * t3250 + t3232 * t3457;
t3116 = -t3148 * t3250 - t3233 * t3451;
t3405 = t3257 * t3519;
t2924 = (-((t3233 * t3405 + t3235 * t3392) * t3544 + ((-t3393 + t3404) * t3251 + pkin(2) * t3405) * t3457 + t2993 * t3235) * t3519 + t3116 * t3245 - (t3233 * t3257 * t3392 + (t3229 * t3235 - t3428 * t3462 - t3235) * t3519) * t3515 + t3567 * t3111) * t3060;
t3558 = 0.2e1 * t2924;
t3113 = t3146 * t3252 + t3232 * t3456;
t3117 = -t3149 * t3252 - t3233 * t3450;
t3403 = t3259 * t3518;
t2925 = (-((t3233 * t3403 + t3235 * t3390) * t3543 + ((-t3391 + t3402) * t3253 + pkin(2) * t3403) * t3456 + t2994 * t3235) * t3518 + t3117 * t3245 - (t3233 * t3259 * t3390 + (t3230 * t3235 - t3426 * t3460 - t3235) * t3518) * t3514 + t3568 * t3113) * t3061;
t3557 = 0.2e1 * t2925;
t3341 = t3066 * t3399;
t3432 = t3245 * t3275;
t3435 = t3235 * t3275;
t3455 = t3233 * t3275;
t3446 = t3235 * t3239;
t3047 = (t3232 * t3446 + t3234 * t3237) * t3542 + t3173 * t3469 + t3576 * t3234;
t3500 = t3047 * t3058;
t3046 = (t3232 * t3237 - t3234 * t3446) * t3542 - t3173 * t3454 + t3576 * t3232;
t3501 = t3046 * t3058;
t3547 = pkin(2) * t3275;
t2953 = t3432 * t3501 - t3349 * t3435 - (-t3235 * t3345 + (-t3138 * t3236 * t3455 + (t3238 * t3547 + t3227) * t3235) * t3517) * t3341 + t3565 * t3275 * t3500;
t3328 = t3005 * t3341;
t3552 = pkin(6) / 0.2e1;
t3556 = -0.2e1 * pkin(2) * t3328 - 0.2e1 * t2953 * t3552;
t3340 = t3082 * t3394;
t3439 = t3235 * t3255;
t3055 = (t3232 * t3439 + t3234 * t3249) * t3541 + t3177 * t3469 + t3575 * t3234;
t3496 = t3055 * t3059;
t3052 = (t3232 * t3249 - t3234 * t3439) * t3541 - t3177 * t3454 + t3575 * t3232;
t3499 = t3052 * t3059;
t2964 = t3432 * t3499 - t3348 * t3435 - (-t3235 * t3344 + (-t3139 * t3248 * t3455 + (t3254 * t3547 + t3228) * t3235) * t3516) * t3340 + t3566 * t3275 * t3496;
t3324 = t3009 * t3340;
t3555 = -0.2e1 * pkin(2) * t3324 - 0.2e1 * t2964 * t3552;
t3339 = t3084 * t3392;
t3438 = t3235 * t3257;
t3056 = (t3232 * t3438 + t3234 * t3251) * t3540 + t3178 * t3469 + t3574 * t3234;
t3495 = t3056 * t3060;
t3053 = (t3232 * t3251 - t3234 * t3438) * t3540 - t3178 * t3454 + t3574 * t3232;
t3498 = t3053 * t3060;
t2965 = t3432 * t3498 - t3347 * t3435 - (-t3235 * t3343 + (-t3140 * t3250 * t3455 + (t3256 * t3547 + t3229) * t3235) * t3515) * t3339 + t3567 * t3275 * t3495;
t3323 = t3010 * t3339;
t3554 = -0.2e1 * pkin(2) * t3323 - 0.2e1 * t2965 * t3552;
t3338 = t3086 * t3390;
t3437 = t3235 * t3259;
t3057 = (t3232 * t3437 + t3234 * t3253) * t3539 + t3179 * t3469 + t3573 * t3234;
t3494 = t3057 * t3061;
t3054 = (t3232 * t3253 - t3234 * t3437) * t3539 - t3179 * t3454 + t3573 * t3232;
t3497 = t3054 * t3061;
t2966 = t3432 * t3497 - t3346 * t3435 - (-t3235 * t3342 + (-t3141 * t3252 * t3455 + (t3258 * t3547 + t3230) * t3235) * t3514) * t3338 + t3568 * t3275 * t3494;
t3322 = t3011 * t3338;
t3553 = -0.2e1 * pkin(2) * t3322 - 0.2e1 * t2966 * t3552;
t3538 = g(3) * t3232;
t3537 = t3246 - g(2);
t3536 = t3247 - g(1);
t3534 = MDP(9) * t3275;
t3533 = t2922 * t3236;
t3532 = t2922 * t3239;
t3531 = t2923 * t3248;
t3530 = t2923 * t3255;
t3529 = t2924 * t3250;
t3528 = t2924 * t3257;
t3527 = t2925 * t3252;
t3526 = t2925 * t3259;
t3525 = t3005 ^ 2 / t3069 ^ 2;
t3523 = t3009 ^ 2 / t3091 ^ 2;
t3522 = t3010 ^ 2 / t3092 ^ 2;
t3521 = t3011 ^ 2 / t3093 ^ 2;
t3070 = t3158 * t3218 - t3317 * t3214;
t3493 = t3058 * t3070;
t3492 = t3058 * t3108;
t3071 = t3159 * t3219 - t3316 * t3215;
t3491 = t3059 * t3071;
t3490 = t3059 * t3115;
t3072 = t3160 * t3220 - t3315 * t3216;
t3489 = t3060 * t3072;
t3488 = t3060 * t3116;
t3073 = t3161 * t3221 - t3314 * t3217;
t3487 = t3061 * t3073;
t3486 = t3061 * t3117;
t3477 = t3293 * t3214;
t3476 = t3293 * t3218;
t3475 = t3109 * t3215;
t3474 = t3109 * t3219;
t3473 = t3111 * t3216;
t3472 = t3111 * t3220;
t3471 = t3113 * t3217;
t3470 = t3113 * t3221;
t3434 = t3236 * t3238;
t3431 = t3248 * t3254;
t3429 = t3250 * t3256;
t3427 = t3252 * t3258;
t3417 = t3058 * t3533;
t3416 = t2922 * t3058 * t3238;
t3415 = t3059 * t3531;
t3414 = t2923 * t3059 * t3254;
t3413 = t3060 * t3529;
t3412 = t2924 * t3060 * t3256;
t3411 = t3061 * t3527;
t3410 = t2925 * t3061 * t3258;
t3276 = 0.1e1 / pkin(3) ^ 2;
t3401 = t3019 ^ 2 / t3096 ^ 2 * t3276;
t3398 = t3023 ^ 2 / t3103 ^ 2 * t3276;
t3397 = t3024 ^ 2 / t3104 ^ 2 * t3276;
t3396 = t3025 ^ 2 / t3105 ^ 2 * t3276;
t3389 = t3047 * t3493;
t3388 = t3214 * t3500;
t3387 = t3218 * t3500;
t3385 = t3055 * t3491;
t3384 = t3215 * t3496;
t3383 = t3219 * t3496;
t3381 = t3056 * t3489;
t3380 = t3216 * t3495;
t3379 = t3220 * t3495;
t3377 = t3057 * t3487;
t3376 = t3217 * t3494;
t3375 = t3221 * t3494;
t3373 = t3293 * t3493;
t3372 = t3058 * t3477;
t3371 = t3058 * t3476;
t3370 = t3109 * t3491;
t3369 = t3059 * t3475;
t3368 = t3059 * t3474;
t3367 = t3111 * t3489;
t3366 = t3060 * t3473;
t3365 = t3060 * t3472;
t3364 = t3113 * t3487;
t3363 = t3061 * t3471;
t3362 = t3061 * t3470;
t3357 = t3047 * t3417;
t3356 = t3047 * t3416;
t3355 = t3055 * t3415;
t3354 = t3055 * t3414;
t3353 = t3056 * t3413;
t3352 = t3056 * t3412;
t3351 = t3057 * t3411;
t3350 = t3057 * t3410;
t3329 = t3058 * t3434 * t3525;
t3327 = t3059 * t3431 * t3523;
t3326 = t3060 * t3429 * t3522;
t3325 = t3061 * t3427 * t3521;
t3321 = -(t3401 + t3525) * t3237 + t3532;
t3320 = -(t3398 + t3523) * t3249 + t3530;
t3319 = -(t3397 + t3522) * t3251 + t3528;
t3318 = -(t3396 + t3521) * t3253 + t3526;
t3313 = t3198 * t3468;
t3312 = t3201 * t3468;
t3311 = t3202 * t3468;
t3310 = t3203 * t3468;
t3309 = 0.2e1 * t3328;
t3308 = 0.2e1 * t3324;
t3307 = 0.2e1 * t3323;
t3306 = 0.2e1 * t3322;
t3297 = t3047 * t3329;
t3296 = t3055 * t3327;
t3295 = t3056 * t3326;
t3294 = t3057 * t3325;
t3292 = pkin(2) * t3560 - pkin(6) * t3401;
t3291 = pkin(2) * t3559 - pkin(6) * t3398;
t3290 = pkin(2) * t3558 - pkin(6) * t3397;
t3289 = pkin(2) * t3557 - pkin(6) * t3396;
t3156 = pkin(2) * t3232 - t3264 * t3454;
t3157 = pkin(2) * t3454 + t3232 * t3264;
t3280 = t3143 * t3546 + (t3156 * t3239 + t3157 * t3237) * t3238;
t3279 = t3147 * t3545 + (t3156 * t3255 + t3157 * t3249) * t3254;
t3278 = t3148 * t3544 + (t3156 * t3257 + t3157 * t3251) * t3256;
t3277 = t3149 * t3543 + (t3156 * t3259 + t3157 * t3253) * t3258;
t3155 = -t3225 * t3244 - t3226 * t3231;
t3154 = -t3225 * t3231 + t3226 * t3244;
t3121 = t3334 * t3234 - t3538;
t3120 = t3335 * t3234 - t3538;
t3119 = t3336 * t3234 - t3538;
t3118 = t3337 * t3234 - t3538;
t3001 = ((-t3161 * t3277 - t3314 * t3592 + (t3161 * t3310 - t3314 * t3577) * t3252) * t3221 + (-t3161 * t3592 + t3314 * t3277 + (-t3161 * t3577 - t3310 * t3314) * t3252) * t3217) * t3086;
t3000 = ((-t3278 * t3160 - t3315 * t3591 + (t3160 * t3311 - t3315 * t3578) * t3250) * t3220 + (-t3160 * t3591 + t3315 * t3278 + (-t3160 * t3578 - t3311 * t3315) * t3250) * t3216) * t3084;
t2999 = ((-t3279 * t3159 - t3316 * t3590 + (t3159 * t3312 - t3316 * t3579) * t3248) * t3219 + (-t3159 * t3590 + t3316 * t3279 + (-t3159 * t3579 - t3312 * t3316) * t3248) * t3215) * t3082;
t2998 = ((-t3280 * t3158 - t3317 * t3589 + (t3158 * t3313 - t3317 * t3580) * t3236) * t3218 + (-t3158 * t3589 + t3280 * t3317 + (-t3158 * t3580 - t3313 * t3317) * t3236) * t3214) * t3066;
t2989 = (-0.2e1 * t3230 + 0.1e1) * t3521;
t2988 = (-0.2e1 * t3229 + 0.1e1) * t3522;
t2987 = (-0.2e1 * t3228 + 0.1e1) * t3523;
t2986 = (-0.2e1 * t3227 + 0.1e1) * t3525;
t2970 = -g(1) * t3217 - g(2) * t3221 + t2973;
t2969 = -g(1) * t3216 - g(2) * t3220 + t2972;
t2968 = -g(1) * t3215 - g(2) * t3219 + t2971;
t2963 = t2970 * t3235 + t3121 * t3233;
t2962 = t2969 * t3235 + t3120 * t3233;
t2961 = t2968 * t3235 + t3119 * t3233;
t2960 = -g(1) * t3214 - g(2) * t3218 + t2967;
t2949 = -t3253 * t3281 + t3585;
t2948 = -t3251 * t3282 + t3586;
t2947 = -t3249 * t3283 + t3587;
t2946 = t2960 * t3235 + t3118 * t3233;
t2942 = -t3237 * t3284 + t3588;
t2941 = t2966 * t3258 - t3252 * t3396;
t2940 = t2966 * t3252 + t3258 * t3396;
t2939 = t2965 * t3256 - t3250 * t3397;
t2938 = t2965 * t3250 + t3256 * t3397;
t2937 = t2964 * t3254 - t3248 * t3398;
t2936 = t2964 * t3248 + t3254 * t3398;
t2935 = t2953 * t3238 - t3236 * t3401;
t2934 = t2953 * t3236 + t3238 * t3401;
t2930 = t2966 * t3253 + t3259 * t3306;
t2929 = t2965 * t3251 + t3257 * t3307;
t2928 = t2964 * t3249 + t3255 * t3308;
t2926 = t2953 * t3237 + t3239 * t3309;
t2921 = t2924 * t3251 + t3257 * t3522;
t2920 = t3251 * t3522 - t3528;
t2919 = t2925 * t3253 + t3259 * t3521;
t2918 = t2923 * t3249 + t3255 * t3523;
t2917 = t3253 * t3521 - t3526;
t2916 = t3249 * t3523 - t3530;
t2915 = t2922 * t3237 + t3239 * t3525;
t2914 = t3237 * t3525 - t3532;
t2913 = (t3258 * t3306 + t3527) * t3252;
t2912 = (t3256 * t3307 + t3529) * t3250;
t2911 = (t3254 * t3308 + t3531) * t3248;
t2910 = (t3238 * t3309 + t3533) * t3236;
t2909 = t3427 * t3557 + (0.4e1 * t3230 - 0.2e1) * t3322;
t2908 = t3429 * t3558 + (0.4e1 * t3229 - 0.2e1) * t3323;
t2907 = t3431 * t3559 + (0.4e1 * t3228 - 0.2e1) * t3324;
t2906 = t3434 * t3560 + (0.4e1 * t3227 - 0.2e1) * t3328;
t2905 = pkin(2) * t3521 - pkin(6) * t2925 + (-t2970 * t3233 + t3121 * t3235) * t3253 + t3585;
t2904 = pkin(2) * t3522 - pkin(6) * t2924 + (-t2969 * t3233 + t3120 * t3235) * t3251 + t3586;
t2903 = pkin(2) * t3523 - pkin(6) * t2923 + (-t2968 * t3233 + t3119 * t3235) * t3249 + t3587;
t2902 = pkin(2) * t3525 - pkin(6) * t2922 + (-t2960 * t3233 + t3118 * t3235) * t3237 + t3588;
t2901 = t2905 * t3258 - t2963 * t3252;
t2900 = t2905 * t3252 + t2963 * t3258;
t2899 = t2904 * t3256 - t2962 * t3250;
t2898 = t2904 * t3250 + t2962 * t3256;
t2897 = t2903 * t3254 - t2961 * t3248;
t2896 = t2903 * t3248 + t2961 * t3254;
t2895 = t2902 * t3238 - t2946 * t3236;
t2894 = t2902 * t3236 + t2946 * t3238;
t2893 = (t2956 + t3289) * t3258 + t3252 * t3553 + t3572 * t3426;
t2892 = (t2955 + t3290) * t3256 + t3250 * t3554 + t3571 * t3428;
t2891 = (t2954 + t3291) * t3254 + t3248 * t3555 + t3570 * t3430;
t2890 = t3258 * t3553 + (-t3289 - t2952) * t3252;
t2889 = t3256 * t3554 + (-t3290 - t2951) * t3250;
t2888 = t3254 * t3555 + (-t3291 - t2950) * t3248;
t2887 = (t2944 + t3292) * t3238 + t3236 * t3556 + t3569 * t3433;
t2886 = t3238 * t3556 + (-t3292 - t2943) * t3236;
t2885 = (-t3252 * t2930 + t3318 * t3258) * t3233 + t3235 * t2941;
t2884 = (-t3258 * t2930 - t3318 * t3252) * t3233 - t3235 * t2940;
t2883 = (-t3250 * t2929 + t3319 * t3256) * t3233 + t3235 * t2939;
t2882 = (-t3256 * t2929 - t3319 * t3250) * t3233 - t3235 * t2938;
t2881 = (-t3248 * t2928 + t3320 * t3254) * t3233 + t3235 * t2937;
t2880 = (-t3254 * t2928 - t3320 * t3248) * t3233 - t3235 * t2936;
t2879 = (-t3236 * t2926 + t3321 * t3238) * t3233 + t3235 * t2935;
t2878 = (-t3238 * t2926 - t3321 * t3236) * t3233 - t3235 * t2934;
t1 = [(t2960 * t3513 + t2968 * t3511 + t2969 * t3510 + t2970 * t3509) * MDP(1) + (-t2922 * t3371 - t2923 * t3368 - t2924 * t3365 - t2925 * t3362) * MDP(2) + (-t2943 * t3371 - t2950 * t3368 - t2951 * t3365 - t2952 * t3362 + (-t2914 * t3513 - t2916 * t3511 - t2917 * t3509 - t2920 * t3510) * t3233) * MDP(3) + (-t2942 * t3371 - t2947 * t3368 - t2948 * t3365 - t2949 * t3362 + (-t2915 * t3513 - t2918 * t3511 - t2919 * t3509 - t2921 * t3510) * t3233) * MDP(4) + (-t2910 * t3371 - t2911 * t3368 - t2912 * t3365 - t2913 * t3362 + (t3218 * t3297 + t3219 * t3296 + t3220 * t3295 + t3221 * t3294) * t3275) * MDP(5) + (-t2906 * t3371 - t2907 * t3368 - t2908 * t3365 - t2909 * t3362 + (-t2986 * t3387 - t2987 * t3383 - t2988 * t3379 - t2989 * t3375) * t3275) * MDP(6) + (-t2934 * t3371 - t2936 * t3368 - t2938 * t3365 - t2940 * t3362 + (-t3218 * t3357 - t3219 * t3355 - t3220 * t3353 - t3221 * t3351) * t3275) * MDP(7) + (-t2935 * t3371 - t2937 * t3368 - t2939 * t3365 - t2941 * t3362 + (-t3218 * t3356 - t3219 * t3354 - t3220 * t3352 - t3221 * t3350) * t3275) * MDP(8) + (-t2953 * t3387 - t2964 * t3383 - t2965 * t3379 - t2966 * t3375) * t3534 + ((t3030 * t2885 - t2893 * t3470) * t3061 + (t2883 * t3029 - t2892 * t3472) * t3060 + (t2881 * t3028 - t2891 * t3474) * t3059 + (t2879 * t3026 - t2887 * t3476) * t3058 + (-t2894 * t3387 - t2896 * t3383 - t2898 * t3379 - t2900 * t3375) * t3275) * MDP(10) + ((t3030 * t2884 - t2890 * t3470) * t3061 + (t2882 * t3029 - t2889 * t3472) * t3060 + (t2880 * t3028 - t2888 * t3474) * t3059 + (t2878 * t3026 - t2886 * t3476) * t3058 + (-t2895 * t3387 - t2897 * t3383 - t2899 * t3379 - t2901 * t3375) * t3275) * MDP(11) + t3155 * MDP(13) - t3154 * MDP(14) + t3536 * MDP(15); (t2960 * t3512 + t2968 * t3508 + t2969 * t3507 + t2970 * t3506) * MDP(1) + (t2922 * t3372 + t2923 * t3369 + t2924 * t3366 + t2925 * t3363) * MDP(2) + (t2943 * t3372 + t2950 * t3369 + t2951 * t3366 + t2952 * t3363 + (-t2914 * t3512 - t2916 * t3508 - t2917 * t3506 - t2920 * t3507) * t3233) * MDP(3) + (t2942 * t3372 + t2947 * t3369 + t2948 * t3366 + t2949 * t3363 + (-t2915 * t3512 - t2918 * t3508 - t2919 * t3506 - t2921 * t3507) * t3233) * MDP(4) + (t2910 * t3372 + t2911 * t3369 + t2912 * t3366 + t2913 * t3363 + (-t3214 * t3297 - t3215 * t3296 - t3216 * t3295 - t3217 * t3294) * t3275) * MDP(5) + (t2906 * t3372 + t2907 * t3369 + t2908 * t3366 + t2909 * t3363 + (t2986 * t3388 + t2987 * t3384 + t2988 * t3380 + t2989 * t3376) * t3275) * MDP(6) + (t2934 * t3372 + t2936 * t3369 + t2938 * t3366 + t2940 * t3363 + (t3214 * t3357 + t3215 * t3355 + t3216 * t3353 + t3217 * t3351) * t3275) * MDP(7) + (t2935 * t3372 + t2937 * t3369 + t2939 * t3366 + t2941 * t3363 + (t3214 * t3356 + t3215 * t3354 + t3216 * t3352 + t3217 * t3350) * t3275) * MDP(8) + (t2953 * t3388 + t2964 * t3384 + t2965 * t3380 + t2966 * t3376) * t3534 + ((t3033 * t2885 + t2893 * t3471) * t3061 + (t3032 * t2883 + t2892 * t3473) * t3060 + (t3031 * t2881 + t2891 * t3475) * t3059 + (t3027 * t2879 + t2887 * t3477) * t3058 + (t2894 * t3388 + t2896 * t3384 + t2898 * t3380 + t2900 * t3376) * t3275) * MDP(10) + ((t3033 * t2884 + t2890 * t3471) * t3061 + (t3032 * t2882 + t2889 * t3473) * t3060 + (t3031 * t2880 + t2888 * t3475) * t3059 + (t3027 * t2878 + t2886 * t3477) * t3058 + (t2895 * t3388 + t2897 * t3384 + t2899 * t3380 + t2901 * t3376) * t3275) * MDP(11) + t3154 * MDP(13) + t3155 * MDP(14) + t3537 * MDP(15); (t2960 * t3505 + t2968 * t3504 + t2969 * t3503 + t2970 * t3502) * MDP(1) + (t2922 * t3492 + t2923 * t3490 + t2924 * t3488 + t2925 * t3486) * MDP(2) + (t2943 * t3492 + t2950 * t3490 + t2951 * t3488 + t2952 * t3486 + (-t2914 * t3505 - t2916 * t3504 - t2917 * t3502 - t2920 * t3503) * t3233) * MDP(3) + (t2942 * t3492 + t2947 * t3490 + t2948 * t3488 + t2949 * t3486 + (-t2915 * t3505 - t2918 * t3504 - t2919 * t3502 - t2921 * t3503) * t3233) * MDP(4) + (t2910 * t3492 + t2911 * t3490 + t2912 * t3488 + t2913 * t3486 + (-t3046 * t3329 - t3052 * t3327 - t3053 * t3326 - t3054 * t3325) * t3275) * MDP(5) + (t2906 * t3492 + t2907 * t3490 + t2908 * t3488 + t2909 * t3486 + (t2986 * t3501 + t2987 * t3499 + t2988 * t3498 + t2989 * t3497) * t3275) * MDP(6) + (t2934 * t3492 + t2936 * t3490 + t2938 * t3488 + t2940 * t3486 + (t3046 * t3417 + t3052 * t3415 + t3053 * t3413 + t3054 * t3411) * t3275) * MDP(7) + (t2935 * t3492 + t2937 * t3490 + t2939 * t3488 + t2941 * t3486 + (t3046 * t3416 + t3052 * t3414 + t3053 * t3412 + t3054 * t3410) * t3275) * MDP(8) + (t2953 * t3501 + t2964 * t3499 + t2965 * t3498 + t2966 * t3497) * t3534 + ((t2885 * t3041 + t2893 * t3117) * t3061 + (t2883 * t3040 + t2892 * t3116) * t3060 + (t2881 * t3039 + t2891 * t3115) * t3059 + (t2879 * t3038 + t2887 * t3108) * t3058 + (t2894 * t3501 + t2896 * t3499 + t2898 * t3498 + t2900 * t3497) * t3275) * MDP(10) + ((t2884 * t3041 + t2890 * t3117) * t3061 + (t2882 * t3040 + t2889 * t3116) * t3060 + (t2880 * t3039 + t2888 * t3115) * t3059 + (t2878 * t3038 + t2886 * t3108) * t3058 + (t2895 * t3501 + t2897 * t3499 + t2899 * t3498 + t2901 * t3497) * t3275) * MDP(11) + (t3245 - g(3)) * MDP(15); (t2960 * t2998 + t2968 * t2999 + t2969 * t3000 + t2970 * t3001) * MDP(1) + (t2922 * t3373 + t2923 * t3370 + t2924 * t3367 + t2925 * t3364) * MDP(2) + (t2943 * t3373 + t2950 * t3370 + t2951 * t3367 + t2952 * t3364 + (-t2914 * t2998 - t2916 * t2999 - t2917 * t3001 - t2920 * t3000) * t3233) * MDP(3) + (t2942 * t3373 + t2947 * t3370 + t2948 * t3367 + t2949 * t3364 + (-t2915 * t2998 - t2918 * t2999 - t2919 * t3001 - t2921 * t3000) * t3233) * MDP(4) + (t2910 * t3373 + t2911 * t3370 + t2912 * t3367 + t2913 * t3364 + (-t3070 * t3297 - t3071 * t3296 - t3072 * t3295 - t3073 * t3294) * t3275) * MDP(5) + (t2906 * t3373 + t2907 * t3370 + t2908 * t3367 + t2909 * t3364 + (t2986 * t3389 + t2987 * t3385 + t2988 * t3381 + t2989 * t3377) * t3275) * MDP(6) + (t2934 * t3373 + t2936 * t3370 + t2938 * t3367 + t2940 * t3364 + (t3070 * t3357 + t3071 * t3355 + t3072 * t3353 + t3073 * t3351) * t3275) * MDP(7) + (t2935 * t3373 + t2937 * t3370 + t2939 * t3367 + t2941 * t3364 + (t3070 * t3356 + t3071 * t3354 + t3072 * t3352 + t3073 * t3350) * t3275) * MDP(8) + (t2953 * t3389 + t2964 * t3385 + t2965 * t3381 + t2966 * t3377) * t3534 + (t2887 * t3373 + t2891 * t3370 + t2892 * t3367 + t2893 * t3364 + t2998 * t2879 + t2999 * t2881 + t3000 * t2883 + t3001 * t2885 + (t2894 * t3389 + t2896 * t3385 + t2898 * t3381 + t2900 * t3377) * t3275) * MDP(10) + (t2886 * t3373 + t2888 * t3370 + t2889 * t3367 + t2890 * t3364 + t2998 * t2878 + t2999 * t2880 + t3000 * t2882 + t3001 * t2884 + (t2895 * t3389 + t2897 * t3385 + t2899 * t3381 + t2901 * t3377) * t3275) * MDP(11) + t3244 * MDP(12) + (-t3536 * t3225 + t3537 * t3226) * MDP(13) + (-t3537 * t3225 - t3536 * t3226) * MDP(14);];
tauX  = t1;
