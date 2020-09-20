% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR2G3P3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [14x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRRRR2G3P3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:13
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR2G3P3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(14,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'P3RRRRR2G3P3A0_gravload_para_pf_mdp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:13:08
% EndTime: 2020-03-09 21:13:10
% DurationCPUTime: 1.66s
% Computational Cost: add. (732->162), mult. (1230->321), div. (303->11), fcn. (1404->30), ass. (0->162)
t3388 = cos(qJ(3,1));
t3471 = t3388 ^ 2;
t3385 = cos(qJ(3,2));
t3470 = t3385 ^ 2;
t3382 = cos(qJ(3,3));
t3469 = t3382 ^ 2;
t3370 = legFrame(3,2);
t3349 = sin(t3370);
t3352 = cos(t3370);
t3336 = t3352 * g(1) - t3349 * g(2);
t3367 = qJ(1,3) + qJ(2,3);
t3343 = sin(t3367);
t3346 = cos(t3367);
t3315 = g(3) * t3346 + t3336 * t3343;
t3374 = sin(qJ(2,3));
t3355 = 0.1e1 / t3374;
t3460 = t3315 * t3355;
t3372 = legFrame(1,2);
t3351 = sin(t3372);
t3354 = cos(t3372);
t3338 = t3354 * g(1) - t3351 * g(2);
t3369 = qJ(1,1) + qJ(2,1);
t3345 = sin(t3369);
t3348 = cos(t3369);
t3319 = g(3) * t3348 + t3338 * t3345;
t3380 = sin(qJ(2,1));
t3357 = 0.1e1 / t3380;
t3459 = t3319 * t3357;
t3371 = legFrame(2,2);
t3350 = sin(t3371);
t3353 = cos(t3371);
t3337 = t3353 * g(1) - t3350 * g(2);
t3368 = qJ(1,2) + qJ(2,2);
t3344 = sin(t3368);
t3347 = cos(t3368);
t3399 = g(3) * t3347 + t3337 * t3344;
t3377 = sin(qJ(2,2));
t3356 = 0.1e1 / t3377;
t3362 = 0.1e1 / t3385;
t3436 = t3356 * t3362;
t3468 = t3399 * t3436;
t3452 = t3344 * t3356;
t3467 = t3399 * t3452;
t3384 = cos(qJ(1,3));
t3466 = pkin(1) * t3384;
t3387 = cos(qJ(1,2));
t3465 = pkin(1) * t3387;
t3390 = cos(qJ(1,1));
t3464 = pkin(1) * t3390;
t3461 = t3399 * t3356;
t3375 = sin(qJ(1,3));
t3383 = cos(qJ(2,3));
t3330 = t3375 * t3374 - t3384 * t3383;
t3458 = t3330 * t3382;
t3378 = sin(qJ(1,2));
t3386 = cos(qJ(2,2));
t3331 = t3378 * t3377 - t3387 * t3386;
t3457 = t3331 * t3385;
t3381 = sin(qJ(1,1));
t3389 = cos(qJ(2,1));
t3332 = t3381 * t3380 - t3390 * t3389;
t3456 = t3332 * t3388;
t3453 = t3343 * t3355;
t3451 = t3345 * t3357;
t3359 = 0.1e1 / t3382;
t3450 = t3349 * t3359;
t3373 = sin(qJ(3,3));
t3449 = t3349 * t3373;
t3448 = t3350 * t3362;
t3376 = sin(qJ(3,2));
t3447 = t3350 * t3376;
t3365 = 0.1e1 / t3388;
t3446 = t3351 * t3365;
t3379 = sin(qJ(3,1));
t3445 = t3351 * t3379;
t3444 = t3352 * t3359;
t3443 = t3352 * t3373;
t3442 = t3353 * t3362;
t3441 = t3353 * t3376;
t3440 = t3354 * t3365;
t3439 = t3354 * t3379;
t3438 = t3355 * t3359;
t3360 = 0.1e1 / t3469;
t3437 = t3355 * t3360;
t3363 = 0.1e1 / t3470;
t3435 = t3356 * t3363;
t3434 = t3357 * t3365;
t3366 = 0.1e1 / t3471;
t3433 = t3357 * t3366;
t3432 = pkin(2) * t3330 * t3469;
t3431 = pkin(2) * t3331 * t3470;
t3430 = pkin(2) * t3332 * t3471;
t3429 = t3383 * t3373 * pkin(1);
t3428 = t3386 * t3376 * pkin(1);
t3427 = t3389 * t3379 * pkin(1);
t3425 = t3376 * t3461;
t3424 = t3373 * t3460;
t3423 = t3379 * t3459;
t3422 = t3315 * t3453;
t3421 = t3315 * t3438;
t3420 = t3315 * t3437;
t3316 = -g(3) * t3343 + t3336 * t3346;
t3419 = t3316 * t3438;
t3418 = t3316 * t3437;
t3416 = t3399 * t3435;
t3318 = -g(3) * t3344 + t3337 * t3347;
t3415 = t3318 * t3436;
t3414 = t3318 * t3435;
t3413 = t3319 * t3451;
t3412 = t3319 * t3434;
t3411 = t3319 * t3433;
t3320 = -g(3) * t3345 + t3338 * t3348;
t3410 = t3320 * t3434;
t3409 = t3320 * t3433;
t3321 = pkin(2) * (t3384 * t3374 + t3375 * t3383) * t3382 + t3375 * pkin(1);
t3408 = t3321 * t3438;
t3322 = pkin(2) * (t3387 * t3377 + t3378 * t3386) * t3385 + t3378 * pkin(1);
t3407 = t3322 * t3436;
t3323 = pkin(2) * (t3390 * t3380 + t3381 * t3389) * t3388 + t3381 * pkin(1);
t3406 = t3323 * t3434;
t3324 = g(3) * t3384 + t3336 * t3375;
t3405 = t3324 * t3438;
t3325 = -g(3) * t3375 + t3336 * t3384;
t3404 = t3325 * t3438;
t3326 = g(3) * t3387 + t3337 * t3378;
t3403 = t3326 * t3436;
t3327 = -g(3) * t3378 + t3337 * t3387;
t3402 = t3327 * t3436;
t3328 = g(3) * t3390 + t3338 * t3381;
t3401 = t3328 * t3434;
t3329 = -g(3) * t3381 + t3338 * t3390;
t3400 = t3329 * t3434;
t3398 = t3362 * t3425;
t3397 = t3363 * t3425;
t3396 = t3359 * t3424;
t3395 = t3360 * t3424;
t3394 = t3365 * t3423;
t3393 = t3366 * t3423;
t3392 = 0.1e1 / pkin(1);
t3391 = 0.1e1 / pkin(2);
t3335 = t3351 * g(1) + t3354 * g(2);
t3334 = t3350 * g(1) + t3353 * g(2);
t3333 = t3349 * g(1) + t3352 * g(2);
t3311 = -t3354 * t3456 + t3445;
t3310 = t3351 * t3456 + t3439;
t3309 = -t3353 * t3457 + t3447;
t3308 = t3350 * t3457 + t3441;
t3307 = -t3352 * t3458 + t3449;
t3306 = t3349 * t3458 + t3443;
t3305 = t3320 * t3388 + t3335 * t3379;
t3304 = t3320 * t3379 - t3335 * t3388;
t3303 = t3318 * t3385 + t3334 * t3376;
t3302 = t3318 * t3376 - t3334 * t3385;
t3301 = t3316 * t3382 + t3333 * t3373;
t3300 = t3316 * t3373 - t3333 * t3382;
t3299 = t3354 * t3430 + (-pkin(2) * t3445 - t3354 * t3464) * t3388 - t3351 * t3427;
t3298 = t3353 * t3431 + (-pkin(2) * t3447 - t3353 * t3465) * t3385 - t3350 * t3428;
t3297 = t3352 * t3432 + (-pkin(2) * t3449 - t3352 * t3466) * t3382 - t3349 * t3429;
t3296 = -t3351 * t3430 + (-pkin(2) * t3439 + t3351 * t3464) * t3388 - t3354 * t3427;
t3295 = -t3350 * t3431 + (-pkin(2) * t3441 + t3350 * t3465) * t3385 - t3353 * t3428;
t3294 = -t3349 * t3432 + (-pkin(2) * t3443 + t3349 * t3466) * t3382 - t3352 * t3429;
t1 = [-g(1) * MDP(14) + ((t3300 * t3450 + t3302 * t3448 + t3304 * t3446) * MDP(12) + (t3301 * t3450 + t3303 * t3448 + t3305 * t3446) * MDP(13)) * t3391 + ((t3307 * t3405 + t3309 * t3403 + t3311 * t3401) * MDP(2) + (t3307 * t3404 + t3309 * t3402 + t3311 * t3400) * MDP(3) + (t3307 * t3421 + t3309 * t3468 + t3311 * t3412) * MDP(5) + (t3307 * t3419 + t3309 * t3415 + t3311 * t3410) * MDP(6) + (t3307 * t3460 + t3309 * t3461 + t3311 * t3459) * MDP(12) + (-t3307 * t3396 - t3309 * t3398 - t3311 * t3394) * MDP(13) + ((t3297 * t3420 + t3298 * t3416 + t3299 * t3411) * MDP(5) + (t3297 * t3418 + t3298 * t3414 + t3299 * t3409) * MDP(6) + (t3297 * t3421 + t3298 * t3468 + t3299 * t3412) * MDP(12) + (-t3297 * t3395 - t3298 * t3397 - t3299 * t3393) * MDP(13)) * t3391) * t3392; -g(2) * MDP(14) + ((t3300 * t3444 + t3302 * t3442 + t3304 * t3440) * MDP(12) + (t3301 * t3444 + t3303 * t3442 + t3305 * t3440) * MDP(13)) * t3391 + ((t3306 * t3405 + t3308 * t3403 + t3310 * t3401) * MDP(2) + (t3306 * t3404 + t3308 * t3402 + t3310 * t3400) * MDP(3) + (t3306 * t3421 + t3308 * t3468 + t3310 * t3412) * MDP(5) + (t3306 * t3419 + t3308 * t3415 + t3310 * t3410) * MDP(6) + (t3306 * t3460 + t3308 * t3461 + t3310 * t3459) * MDP(12) + (-t3306 * t3396 - t3308 * t3398 - t3310 * t3394) * MDP(13) + ((t3294 * t3420 + t3295 * t3416 + t3296 * t3411) * MDP(5) + (t3294 * t3418 + t3295 * t3414 + t3296 * t3409) * MDP(6) + (t3294 * t3421 + t3295 * t3468 + t3296 * t3412) * MDP(12) + (-t3294 * t3395 - t3295 * t3397 - t3296 * t3393) * MDP(13)) * t3391) * t3392; -g(3) * MDP(14) + ((-t3324 * t3453 - t3326 * t3452 - t3328 * t3451) * MDP(2) + (-t3325 * t3453 - t3327 * t3452 - t3329 * t3451) * MDP(3) + (-t3413 - t3422 - t3467) * MDP(5) + (-t3316 * t3453 - t3318 * t3452 - t3320 * t3451) * MDP(6) + (-t3382 * t3422 - t3385 * t3467 - t3388 * t3413) * MDP(12) + (t3343 * t3424 + t3344 * t3425 + t3345 * t3423) * MDP(13) + ((t3315 * t3408 + t3319 * t3406 + t3399 * t3407) * MDP(5) + (t3316 * t3408 + t3318 * t3407 + t3320 * t3406) * MDP(6) + (t3321 * t3460 + t3322 * t3461 + t3323 * t3459) * MDP(12) + (-t3321 * t3396 - t3322 * t3398 - t3323 * t3394) * MDP(13)) * t3391) * t3392;];
taugX  = t1;
