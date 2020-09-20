% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR2G2A0
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
%   see P3RRRRR2G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:10
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR2G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(14,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G2A0_gravload_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'P3RRRRR2G2A0_gravload_para_pf_mdp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:09:59
% EndTime: 2020-03-09 21:10:01
% DurationCPUTime: 1.64s
% Computational Cost: add. (810->177), mult. (1215->316), div. (294->14), fcn. (1368->42), ass. (0->169)
t3404 = cos(qJ(3,1));
t3499 = t3404 ^ 2;
t3401 = cos(qJ(3,2));
t3498 = t3401 ^ 2;
t3398 = cos(qJ(3,3));
t3497 = t3398 ^ 2;
t3388 = legFrame(1,2);
t3367 = sin(t3388);
t3370 = cos(t3388);
t3348 = -t3370 * g(1) + t3367 * g(2);
t3385 = qJ(1,1) + qJ(2,1);
t3361 = sin(t3385);
t3364 = cos(t3385);
t3415 = g(3) * t3361 + t3348 * t3364;
t3396 = sin(qJ(2,1));
t3373 = 0.1e1 / t3396;
t3381 = 0.1e1 / t3404;
t3449 = t3373 * t3381;
t3496 = t3415 * t3449;
t3387 = legFrame(2,2);
t3366 = sin(t3387);
t3369 = cos(t3387);
t3346 = -t3369 * g(1) + t3366 * g(2);
t3384 = qJ(1,2) + qJ(2,2);
t3360 = sin(t3384);
t3363 = cos(t3384);
t3416 = g(3) * t3360 + t3346 * t3363;
t3393 = sin(qJ(2,2));
t3372 = 0.1e1 / t3393;
t3378 = 0.1e1 / t3401;
t3451 = t3372 * t3378;
t3495 = t3416 * t3451;
t3386 = legFrame(3,2);
t3365 = sin(t3386);
t3368 = cos(t3386);
t3344 = -t3368 * g(1) + t3365 * g(2);
t3383 = qJ(1,3) + qJ(2,3);
t3359 = sin(t3383);
t3362 = cos(t3383);
t3417 = g(3) * t3359 + t3344 * t3362;
t3390 = sin(qJ(2,3));
t3371 = 0.1e1 / t3390;
t3375 = 0.1e1 / t3398;
t3453 = t3371 * t3375;
t3494 = t3417 * t3453;
t3493 = -2 * pkin(1);
t3391 = sin(qJ(1,3));
t3492 = pkin(1) * t3391;
t3394 = sin(qJ(1,2));
t3491 = pkin(1) * t3394;
t3397 = sin(qJ(1,1));
t3490 = pkin(1) * t3397;
t3489 = qJ(2,1) - qJ(3,1);
t3488 = qJ(2,1) + qJ(3,1);
t3487 = qJ(2,2) - qJ(3,2);
t3486 = qJ(2,2) + qJ(3,2);
t3485 = qJ(2,3) - qJ(3,3);
t3484 = qJ(2,3) + qJ(3,3);
t3483 = t3417 * t3371;
t3482 = t3417 * t3398;
t3481 = t3416 * t3372;
t3480 = t3416 * t3401;
t3479 = t3415 * t3373;
t3478 = t3415 * t3404;
t3389 = sin(qJ(3,3));
t3477 = t3417 * t3389;
t3392 = sin(qJ(3,2));
t3476 = t3416 * t3392;
t3395 = sin(qJ(3,1));
t3475 = t3415 * t3395;
t3400 = cos(qJ(1,3));
t3474 = (t3400 * t3493 + (-cos(qJ(1,3) + t3485) - cos(qJ(1,3) + t3484)) * pkin(2)) / (sin(t3484) + sin(t3485));
t3403 = cos(qJ(1,2));
t3473 = (t3403 * t3493 + (-cos(qJ(1,2) + t3487) - cos(qJ(1,2) + t3486)) * pkin(2)) / (sin(t3486) + sin(t3487));
t3406 = cos(qJ(1,1));
t3472 = (t3406 * t3493 + (-cos(qJ(1,1) + t3489) - cos(qJ(1,1) + t3488)) * pkin(2)) / (sin(t3488) + sin(t3489));
t3399 = cos(qJ(2,3));
t3341 = t3400 * t3390 + t3391 * t3399;
t3471 = t3341 * t3398;
t3402 = cos(qJ(2,2));
t3342 = t3403 * t3393 + t3394 * t3402;
t3470 = t3342 * t3401;
t3405 = cos(qJ(2,1));
t3343 = t3406 * t3396 + t3397 * t3405;
t3469 = t3343 * t3404;
t3468 = t3362 * t3371;
t3467 = t3363 * t3372;
t3466 = t3364 * t3373;
t3465 = t3365 * t3375;
t3464 = t3365 * t3389;
t3463 = t3366 * t3378;
t3462 = t3366 * t3392;
t3461 = t3367 * t3381;
t3460 = t3367 * t3395;
t3459 = t3368 * t3375;
t3458 = t3368 * t3389;
t3457 = t3369 * t3378;
t3456 = t3369 * t3392;
t3455 = t3370 * t3381;
t3454 = t3370 * t3395;
t3376 = 0.1e1 / t3497;
t3452 = t3371 * t3376;
t3379 = 0.1e1 / t3498;
t3450 = t3372 * t3379;
t3382 = 0.1e1 / t3499;
t3448 = t3373 * t3382;
t3447 = t3341 * t3497 * pkin(2);
t3446 = t3342 * t3498 * pkin(2);
t3445 = t3343 * t3499 * pkin(2);
t3444 = t3399 * t3389 * pkin(1);
t3443 = t3402 * t3392 * pkin(1);
t3442 = t3405 * t3395 * pkin(1);
t3438 = t3371 * t3477;
t3437 = t3372 * t3476;
t3436 = t3373 * t3475;
t3434 = t3417 * t3452;
t3432 = t3416 * t3450;
t3430 = t3415 * t3448;
t3329 = g(3) * t3362 - t3344 * t3359;
t3429 = t3329 * t3453;
t3428 = t3329 * t3452;
t3330 = g(3) * t3363 - t3346 * t3360;
t3427 = t3330 * t3451;
t3426 = t3330 * t3450;
t3331 = g(3) * t3364 - t3348 * t3361;
t3425 = t3331 * t3449;
t3424 = t3331 * t3448;
t3332 = g(3) * t3391 + t3344 * t3400;
t3423 = t3332 * t3453;
t3333 = g(3) * t3394 + t3346 * t3403;
t3422 = t3333 * t3451;
t3334 = g(3) * t3397 + t3348 * t3406;
t3421 = t3334 * t3449;
t3335 = g(3) * t3400 - t3344 * t3391;
t3420 = t3335 * t3453;
t3336 = g(3) * t3403 - t3346 * t3394;
t3419 = t3336 * t3451;
t3337 = g(3) * t3406 - t3348 * t3397;
t3418 = t3337 * t3449;
t3414 = t3375 * t3438;
t3413 = t3376 * t3438;
t3412 = t3378 * t3437;
t3411 = t3379 * t3437;
t3410 = t3381 * t3436;
t3409 = t3382 * t3436;
t3408 = 1 / pkin(1);
t3407 = 0.1e1 / pkin(2);
t3349 = t3367 * g(1) + t3370 * g(2);
t3347 = t3366 * g(1) + t3369 * g(2);
t3345 = t3365 * g(1) + t3368 * g(2);
t3319 = t3370 * t3469 + t3460;
t3318 = -t3367 * t3469 + t3454;
t3317 = t3369 * t3470 + t3462;
t3316 = -t3366 * t3470 + t3456;
t3315 = t3368 * t3471 + t3464;
t3314 = -t3365 * t3471 + t3458;
t3313 = t3331 * t3404 + t3349 * t3395;
t3312 = t3331 * t3395 - t3349 * t3404;
t3311 = t3330 * t3401 + t3347 * t3392;
t3310 = t3330 * t3392 - t3347 * t3401;
t3309 = t3329 * t3398 + t3345 * t3389;
t3308 = t3329 * t3389 - t3345 * t3398;
t3307 = -t3370 * t3445 + (-pkin(2) * t3460 - t3370 * t3490) * t3404 - t3367 * t3442;
t3306 = t3367 * t3445 + (-pkin(2) * t3454 + t3367 * t3490) * t3404 - t3370 * t3442;
t3305 = -t3369 * t3446 + (-pkin(2) * t3462 - t3369 * t3491) * t3401 - t3366 * t3443;
t3304 = t3366 * t3446 + (-pkin(2) * t3456 + t3366 * t3491) * t3401 - t3369 * t3443;
t3303 = -t3368 * t3447 + (-pkin(2) * t3464 - t3368 * t3492) * t3398 - t3365 * t3444;
t3302 = t3365 * t3447 + (-pkin(2) * t3458 + t3365 * t3492) * t3398 - t3368 * t3444;
t1 = [-g(1) * MDP(14) + ((t3308 * t3465 + t3310 * t3463 + t3312 * t3461) * MDP(12) + (t3309 * t3465 + t3311 * t3463 + t3313 * t3461) * MDP(13)) * t3407 + ((t3315 * t3423 + t3317 * t3422 + t3319 * t3421) * MDP(2) + (t3315 * t3420 + t3317 * t3419 + t3319 * t3418) * MDP(3) + (t3315 * t3494 + t3317 * t3495 + t3319 * t3496) * MDP(5) + (t3315 * t3429 + t3317 * t3427 + t3319 * t3425) * MDP(6) + (t3315 * t3483 + t3317 * t3481 + t3319 * t3479) * MDP(12) + (-t3315 * t3414 - t3317 * t3412 - t3319 * t3410) * MDP(13) + ((t3303 * t3434 + t3305 * t3432 + t3307 * t3430) * MDP(5) + (t3303 * t3428 + t3305 * t3426 + t3307 * t3424) * MDP(6) + (t3303 * t3494 + t3305 * t3495 + t3307 * t3496) * MDP(12) + (-t3303 * t3413 - t3305 * t3411 - t3307 * t3409) * MDP(13)) * t3407) * t3408; -g(2) * MDP(14) + ((t3308 * t3459 + t3310 * t3457 + t3312 * t3455) * MDP(12) + (t3309 * t3459 + t3311 * t3457 + t3313 * t3455) * MDP(13)) * t3407 + ((t3314 * t3423 + t3316 * t3422 + t3318 * t3421) * MDP(2) + (t3314 * t3420 + t3316 * t3419 + t3318 * t3418) * MDP(3) + (t3314 * t3494 + t3316 * t3495 + t3318 * t3496) * MDP(5) + (t3314 * t3429 + t3316 * t3427 + t3318 * t3425) * MDP(6) + (t3314 * t3483 + t3316 * t3481 + t3318 * t3479) * MDP(12) + (-t3314 * t3414 - t3316 * t3412 - t3318 * t3410) * MDP(13) + ((t3302 * t3434 + t3304 * t3432 + t3306 * t3430) * MDP(5) + (t3302 * t3428 + t3304 * t3426 + t3306 * t3424) * MDP(6) + (t3302 * t3494 + t3304 * t3495 + t3306 * t3496) * MDP(12) + (-t3302 * t3413 - t3304 * t3411 - t3306 * t3409) * MDP(13)) * t3407) * t3408; -g(3) * MDP(14) + ((t3332 * t3468 + t3333 * t3467 + t3334 * t3466) * MDP(2) + (t3335 * t3468 + t3336 * t3467 + t3337 * t3466) * MDP(3) + (t3415 * t3466 + t3416 * t3467 + t3417 * t3468) * MDP(5) + (t3329 * t3468 + t3330 * t3467 + t3331 * t3466) * MDP(6) + (t3466 * t3478 + t3467 * t3480 + t3468 * t3482) * MDP(12) + (-t3362 * t3438 - t3363 * t3437 - t3364 * t3436) * MDP(13) + ((t3415 * t3472 + t3416 * t3473 + t3417 * t3474) * MDP(5) + (t3329 * t3474 + t3330 * t3473 + t3331 * t3472) * MDP(6) + (t3472 * t3478 + t3473 * t3480 + t3474 * t3482) * MDP(12) + (-t3472 * t3475 - t3473 * t3476 - t3474 * t3477) * MDP(13)) * t3407) * t3408;];
taugX  = t1;
