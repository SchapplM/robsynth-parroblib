% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4PRRRR1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P4PRRRR1G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [4x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR1G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_gravload_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_gravload_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:59:59
% EndTime: 2020-08-07 11:00:01
% DurationCPUTime: 1.70s
% Computational Cost: add. (453->137), mult. (1077->277), div. (272->13), fcn. (1176->26), ass. (0->135)
t3293 = legFrame(1,2);
t3260 = sin(t3293);
t3264 = cos(t3293);
t3253 = t3264 * g(1) - t3260 * g(2);
t3298 = sin(qJ(3,1));
t3304 = cos(qJ(3,1));
t3249 = t3260 * g(1) + t3264 * g(2);
t3299 = sin(qJ(2,1));
t3305 = cos(qJ(2,1));
t3340 = g(3) * t3305 + t3249 * t3299;
t3220 = t3253 * t3304 + t3298 * t3340;
t3221 = -t3253 * t3298 + t3304 * t3340;
t3229 = -g(3) * t3299 + t3249 * t3305;
t3284 = 0.1e1 / t3304;
t3285 = 0.1e1 / t3304 ^ 2;
t3352 = t3229 * t3284 * t3298;
t3279 = 0.1e1 / t3299;
t3364 = t3279 * t3305;
t3376 = t3229 * t3298 ^ 2;
t3399 = (-MDP(10) * t3352 + (t3298 * (-MDP(3) * t3229 + MDP(4) * t3340) + MDP(11) * t3376) * t3285) * t3364 - (MDP(10) * t3220 + t3221 * MDP(11)) * t3284;
t3292 = legFrame(2,2);
t3259 = sin(t3292);
t3263 = cos(t3292);
t3252 = t3263 * g(1) - t3259 * g(2);
t3296 = sin(qJ(3,2));
t3302 = cos(qJ(3,2));
t3248 = t3259 * g(1) + t3263 * g(2);
t3297 = sin(qJ(2,2));
t3303 = cos(qJ(2,2));
t3341 = g(3) * t3303 + t3248 * t3297;
t3218 = t3252 * t3302 + t3296 * t3341;
t3219 = -t3252 * t3296 + t3302 * t3341;
t3228 = -g(3) * t3297 + t3248 * t3303;
t3282 = 0.1e1 / t3302;
t3283 = 0.1e1 / t3302 ^ 2;
t3353 = t3228 * t3282 * t3296;
t3277 = 0.1e1 / t3297;
t3366 = t3277 * t3303;
t3377 = t3228 * t3296 ^ 2;
t3398 = (-MDP(10) * t3353 + (t3296 * (-MDP(3) * t3228 + MDP(4) * t3341) + MDP(11) * t3377) * t3283) * t3366 - (MDP(10) * t3218 + t3219 * MDP(11)) * t3282;
t3291 = legFrame(3,2);
t3258 = sin(t3291);
t3262 = cos(t3291);
t3251 = t3262 * g(1) - t3258 * g(2);
t3294 = sin(qJ(3,3));
t3300 = cos(qJ(3,3));
t3247 = t3258 * g(1) + t3262 * g(2);
t3295 = sin(qJ(2,3));
t3301 = cos(qJ(2,3));
t3342 = g(3) * t3301 + t3247 * t3295;
t3216 = t3251 * t3300 + t3294 * t3342;
t3217 = -t3251 * t3294 + t3300 * t3342;
t3227 = -g(3) * t3295 + t3247 * t3301;
t3280 = 0.1e1 / t3300;
t3281 = 0.1e1 / t3300 ^ 2;
t3354 = t3227 * t3280 * t3294;
t3275 = 0.1e1 / t3295;
t3368 = t3275 * t3301;
t3378 = t3227 * t3294 ^ 2;
t3397 = (-MDP(10) * t3354 + (t3294 * (-MDP(3) * t3227 + MDP(4) * t3342) + MDP(11) * t3378) * t3281) * t3368 - (MDP(10) * t3216 + t3217 * MDP(11)) * t3280;
t3290 = legFrame(4,2);
t3257 = sin(t3290);
t3261 = cos(t3290);
t3250 = t3261 * g(1) - t3257 * g(2);
t3286 = sin(qJ(3,4));
t3288 = cos(qJ(3,4));
t3246 = t3257 * g(1) + t3261 * g(2);
t3287 = sin(qJ(2,4));
t3289 = cos(qJ(2,4));
t3343 = g(3) * t3289 + t3246 * t3287;
t3214 = t3250 * t3288 + t3286 * t3343;
t3215 = -t3250 * t3286 + t3288 * t3343;
t3223 = -g(3) * t3287 + t3246 * t3289;
t3272 = 0.1e1 / t3288;
t3273 = 0.1e1 / t3288 ^ 2;
t3355 = t3223 * t3272 * t3286;
t3271 = 0.1e1 / t3287;
t3370 = t3271 * t3289;
t3379 = t3223 * t3286 ^ 2;
t3396 = (-MDP(10) * t3355 + (t3286 * (-MDP(3) * t3223 + MDP(4) * t3343) + MDP(11) * t3379) * t3273) * t3370 - (MDP(10) * t3214 + t3215 * MDP(11)) * t3272;
t3371 = t3271 * t3272;
t3369 = t3275 * t3280;
t3367 = t3277 * t3282;
t3365 = t3279 * t3284;
t3363 = t3287 * t3288;
t3362 = t3295 * t3300;
t3361 = t3297 * t3302;
t3360 = t3299 * t3304;
t3351 = t3246 * t3371;
t3350 = t3247 * t3369;
t3349 = t3248 * t3367;
t3348 = t3249 * t3365;
t3347 = t3286 * t3370;
t3346 = t3294 * t3368;
t3345 = t3296 * t3366;
t3344 = t3298 * t3364;
t3306 = xP(4);
t3268 = sin(t3306);
t3269 = cos(t3306);
t3310 = koppelP(1,2);
t3314 = koppelP(1,1);
t3316 = -t3268 * t3310 + t3269 * t3314;
t3317 = t3268 * t3314 + t3269 * t3310;
t3210 = t3260 * t3316 + t3264 * t3317;
t3331 = t3210 * t3285 * t3344;
t3307 = koppelP(4,2);
t3311 = koppelP(4,1);
t3322 = -t3268 * t3307 + t3269 * t3311;
t3323 = t3268 * t3311 + t3269 * t3307;
t3211 = t3257 * t3322 + t3261 * t3323;
t3330 = t3211 * t3273 * t3347;
t3308 = koppelP(3,2);
t3312 = koppelP(3,1);
t3320 = -t3268 * t3308 + t3269 * t3312;
t3321 = t3268 * t3312 + t3269 * t3308;
t3212 = t3258 * t3320 + t3262 * t3321;
t3329 = t3212 * t3281 * t3346;
t3309 = koppelP(2,2);
t3313 = koppelP(2,1);
t3318 = -t3268 * t3309 + t3269 * t3313;
t3319 = t3268 * t3313 + t3269 * t3309;
t3213 = t3259 * t3318 + t3263 * t3319;
t3328 = t3213 * t3283 * t3345;
t3315 = 0.1e1 / pkin(2);
t3255 = t3269 * g(1) + t3268 * g(2);
t3254 = t3268 * g(1) - t3269 * g(2);
t3245 = t3260 * t3298 + t3264 * t3360;
t3244 = t3259 * t3296 + t3263 * t3361;
t3243 = t3258 * t3294 + t3262 * t3362;
t3242 = t3260 * t3360 - t3264 * t3298;
t3241 = t3259 * t3361 - t3263 * t3296;
t3240 = t3258 * t3362 - t3262 * t3294;
t3239 = t3257 * t3286 + t3261 * t3363;
t3238 = t3257 * t3363 - t3261 * t3286;
t1 = [(-t3238 * t3351 - t3240 * t3350 - t3241 * t3349 - t3242 * t3348) * MDP(1) + (-t3268 * t3254 - t3269 * t3255) * MDP(15) + (t3396 * t3261 + t3397 * t3262 + t3398 * t3263 + t3399 * t3264) * t3315; (-t3239 * t3351 - t3243 * t3350 - t3244 * t3349 - t3245 * t3348) * MDP(1) + (t3269 * t3254 - t3268 * t3255) * MDP(15) + (-t3396 * t3257 - t3397 * t3258 - t3398 * t3259 - t3399 * t3260) * t3315; (-t3246 * t3370 - t3247 * t3368 - t3248 * t3366 - t3249 * t3364) * MDP(1) - g(3) * MDP(15) + ((t3223 * t3371 + t3227 * t3369 + t3228 * t3367 + t3229 * t3365) * MDP(3) + (-t3340 * t3365 - t3341 * t3367 - t3342 * t3369 - t3343 * t3371) * MDP(4) + (t3223 * t3271 + t3227 * t3275 + t3228 * t3277 + t3229 * t3279) * MDP(10) + (-t3271 * t3355 - t3275 * t3354 - t3277 * t3353 - t3279 * t3352) * MDP(11)) * t3315; (-(-t3242 * t3317 + t3245 * t3316) * t3348 - (-t3241 * t3319 + t3244 * t3318) * t3349 - (-t3240 * t3321 + t3243 * t3320) * t3350 - (-t3238 * t3323 + t3239 * t3322) * t3351) * MDP(1) + t3254 * MDP(13) + t3255 * MDP(14) + ((t3223 * t3330 + t3227 * t3329 + t3228 * t3328 + t3229 * t3331) * MDP(3) + (-t3328 * t3341 - t3329 * t3342 - t3330 * t3343 - t3331 * t3340) * MDP(4) + ((t3228 * t3345 + t3218) * t3282 * t3213 + (t3227 * t3346 + t3216) * t3280 * t3212 + (t3223 * t3347 + t3214) * t3272 * t3211 + (t3229 * t3344 + t3220) * t3284 * t3210) * MDP(10) + ((-t3283 * t3366 * t3377 + t3282 * t3219) * t3213 + (-t3281 * t3368 * t3378 + t3280 * t3217) * t3212 + (-t3273 * t3370 * t3379 + t3272 * t3215) * t3211 + (-t3285 * t3364 * t3376 + t3284 * t3221) * t3210) * MDP(11)) * t3315;];
taugX  = t1;
