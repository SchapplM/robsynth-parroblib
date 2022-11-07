% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR8V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-04 17:05
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V1G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:05:38
% EndTime: 2022-11-04 17:05:40
% DurationCPUTime: 1.85s
% Computational Cost: add. (984->181), mult. (1776->305), div. (159->9), fcn. (1668->26), ass. (0->148)
t3448 = sin(qJ(1,3));
t3454 = cos(qJ(1,3));
t3444 = legFrame(3,2);
t3428 = sin(t3444);
t3431 = cos(t3444);
t3472 = t3431 * g(1) - t3428 * g(2);
t3395 = g(3) * t3448 - t3472 * t3454;
t3441 = pkin(4) + qJ(3,3);
t3437 = 0.1e1 / t3441;
t3510 = t3437 * t3454;
t3496 = t3395 * t3510;
t3421 = cos(pkin(5)) * pkin(2) + pkin(1);
t3447 = sin(qJ(2,3));
t3453 = cos(qJ(2,3));
t3539 = pkin(2) * sin(pkin(5));
t3469 = t3421 * t3453 - t3447 * t3539;
t3525 = 0.1e1 / t3469 * t3437;
t3552 = t3395 * t3525;
t3450 = sin(qJ(1,2));
t3456 = cos(qJ(1,2));
t3445 = legFrame(2,2);
t3429 = sin(t3445);
t3432 = cos(t3445);
t3471 = t3432 * g(1) - t3429 * g(2);
t3396 = g(3) * t3450 - t3471 * t3456;
t3442 = pkin(4) + qJ(3,2);
t3438 = 0.1e1 / t3442;
t3509 = t3438 * t3456;
t3551 = t3396 * t3509;
t3449 = sin(qJ(2,2));
t3455 = cos(qJ(2,2));
t3468 = t3421 * t3455 - t3449 * t3539;
t3524 = 0.1e1 / t3468 * t3438;
t3550 = t3396 * t3524;
t3443 = pkin(4) + qJ(3,1);
t3439 = 0.1e1 / t3443;
t3458 = cos(qJ(1,1));
t3508 = t3439 * t3458;
t3452 = sin(qJ(1,1));
t3446 = legFrame(1,2);
t3430 = sin(t3446);
t3433 = cos(t3446);
t3470 = t3433 * g(1) - t3430 * g(2);
t3544 = g(3) * t3452 - t3470 * t3458;
t3549 = t3544 * t3508;
t3451 = sin(qJ(2,1));
t3457 = cos(qJ(2,1));
t3467 = t3421 * t3457 - t3451 * t3539;
t3523 = 0.1e1 / t3467 * t3439;
t3548 = t3544 * t3523;
t3547 = MDP(3) - MDP(13);
t3543 = MDP(14) * pkin(1);
t3434 = qJ(2,3) + pkin(5);
t3425 = cos(t3434);
t3542 = pkin(2) * t3425;
t3435 = qJ(2,2) + pkin(5);
t3426 = cos(t3435);
t3541 = pkin(2) * t3426;
t3436 = qJ(2,1) + pkin(5);
t3427 = cos(t3436);
t3540 = pkin(2) * t3427;
t3531 = t3453 * pkin(1);
t3530 = t3455 * pkin(1);
t3529 = t3457 * pkin(1);
t3528 = t3395 * t3437;
t3527 = t3396 * t3438;
t3526 = t3544 * t3439;
t3413 = 0.1e1 / (t3531 + t3542);
t3522 = t3413 * t3428;
t3521 = t3413 * t3431;
t3414 = 0.1e1 / (t3530 + t3541);
t3520 = t3414 * t3429;
t3519 = t3414 * t3432;
t3415 = 0.1e1 / (t3529 + t3540);
t3518 = t3415 * t3430;
t3517 = t3415 * t3433;
t3516 = t3421 * t3431;
t3515 = t3421 * t3432;
t3514 = t3421 * t3433;
t3513 = t3428 * t3421;
t3512 = t3429 * t3421;
t3511 = t3430 * t3421;
t3507 = t3431 * t3539;
t3506 = t3432 * t3539;
t3505 = t3433 * t3539;
t3504 = t3428 * t3539;
t3503 = t3429 * t3539;
t3502 = t3430 * t3539;
t3466 = g(3) * t3454 + t3472 * t3448;
t3492 = t3466 * t3525;
t3465 = g(3) * t3456 + t3471 * t3450;
t3491 = t3465 * t3524;
t3464 = g(3) * t3458 + t3470 * t3452;
t3490 = t3464 * t3523;
t3418 = t3454 * t3531;
t3384 = g(3) * (-t3454 * qJ(3,3) + t3448 * t3531) - t3472 * (t3448 * qJ(3,3) + t3418);
t3489 = t3384 * t3525;
t3419 = t3456 * t3530;
t3385 = g(3) * (-t3456 * qJ(3,2) + t3450 * t3530) - t3471 * (t3450 * qJ(3,2) + t3419);
t3488 = t3385 * t3524;
t3420 = t3458 * t3529;
t3383 = g(3) * (-t3458 * qJ(3,1) + t3452 * t3529) - t3470 * (t3452 * qJ(3,1) + t3420);
t3487 = t3383 * t3523;
t3484 = t3427 * t3548;
t3483 = t3457 * t3548;
t3482 = t3426 * t3550;
t3481 = t3449 * t3550;
t3480 = t3455 * t3550;
t3422 = sin(t3434);
t3479 = t3422 * t3552;
t3478 = t3425 * t3552;
t3477 = t3447 * t3552;
t3476 = t3453 * t3552;
t3423 = sin(t3435);
t3475 = t3423 * t3550;
t3424 = sin(t3436);
t3474 = t3424 * t3548;
t3473 = t3451 * t3548;
t3410 = t3428 * g(1) + t3431 * g(2);
t3377 = -t3410 * t3453 + t3447 * t3466;
t3411 = t3429 * g(1) + t3432 * g(2);
t3379 = -t3411 * t3455 + t3449 * t3465;
t3412 = t3430 * g(1) + t3433 * g(2);
t3381 = -t3412 * t3457 + t3451 * t3464;
t3463 = t3377 * t3522 + t3379 * t3520 + t3381 * t3518;
t3462 = t3377 * t3521 + t3379 * t3519 + t3381 * t3517;
t3406 = t3451 * t3421 + t3457 * t3539;
t3405 = t3449 * t3421 + t3455 * t3539;
t3404 = t3447 * t3421 + t3453 * t3539;
t3388 = -t3458 * t3443 + t3467 * t3452;
t3387 = -t3456 * t3442 + t3468 * t3450;
t3386 = -t3454 * t3441 + t3469 * t3448;
t3382 = t3412 * t3451 + t3457 * t3464;
t3380 = t3411 * t3449 + t3455 * t3465;
t3378 = t3410 * t3447 + t3453 * t3466;
t3376 = t3412 * t3424 + t3427 * t3464;
t3375 = -t3412 * t3427 + t3424 * t3464;
t3374 = t3411 * t3423 + t3426 * t3465;
t3373 = -t3411 * t3426 + t3423 * t3465;
t3372 = t3410 * t3422 + t3425 * t3466;
t3371 = -t3410 * t3425 + t3422 * t3466;
t3370 = (t3452 * t3514 + t3502) * t3457 + (-t3452 * t3505 + t3511) * t3451;
t3369 = (t3450 * t3515 + t3503) * t3455 + (-t3450 * t3506 + t3512) * t3449;
t3368 = (t3448 * t3516 + t3504) * t3453 + (-t3448 * t3507 + t3513) * t3447;
t3367 = (-t3452 * t3511 + t3505) * t3457 + (t3452 * t3502 + t3514) * t3451;
t3366 = (-t3450 * t3512 + t3506) * t3455 + (t3450 * t3503 + t3515) * t3449;
t3365 = (-t3448 * t3513 + t3507) * t3453 + (t3448 * t3504 + t3516) * t3447;
t1 = [(t3368 * t3552 + t3369 * t3550 + t3370 * t3548) * MDP(2) + (t3368 * t3476 + t3369 * t3480 + t3370 * t3483 + t3463) * MDP(9) + (-t3368 * t3477 - t3369 * t3481 - t3370 * t3473 + t3378 * t3522 + t3380 * t3520 + t3382 * t3518) * MDP(10) + (t3368 * t3478 + t3369 * t3482 + t3370 * t3484 + t3371 * t3522 + t3373 * t3520 + t3375 * t3518) * MDP(11) + (-t3368 * t3479 - t3369 * t3475 - t3370 * t3474 + t3372 * t3522 + t3374 * t3520 + t3376 * t3518) * MDP(12) + (t3370 * t3487 - (t3388 * t3433 + t3430 * t3406) * t3526 + t3369 * t3488 - (t3387 * t3432 + t3429 * t3405) * t3527 + t3368 * t3489 - (t3386 * t3431 + t3428 * t3404) * t3528) * MDP(14) - g(1) * MDP(15) + t3463 * t3543 + t3547 * (t3368 * t3492 + t3369 * t3491 + t3370 * t3490); (t3365 * t3552 + t3366 * t3550 + t3367 * t3548) * MDP(2) + (t3365 * t3476 + t3366 * t3480 + t3367 * t3483 + t3462) * MDP(9) + (-t3365 * t3477 - t3366 * t3481 - t3367 * t3473 + t3378 * t3521 + t3380 * t3519 + t3382 * t3517) * MDP(10) + (t3365 * t3478 + t3366 * t3482 + t3367 * t3484 + t3371 * t3521 + t3373 * t3519 + t3375 * t3517) * MDP(11) + (-t3365 * t3479 - t3366 * t3475 - t3367 * t3474 + t3372 * t3521 + t3374 * t3519 + t3376 * t3517) * MDP(12) + (t3367 * t3487 - (-t3388 * t3430 + t3433 * t3406) * t3526 + t3366 * t3488 - (-t3387 * t3429 + t3432 * t3405) * t3527 + t3365 * t3489 - (-t3386 * t3428 + t3431 * t3404) * t3528) * MDP(14) - g(2) * MDP(15) + t3462 * t3543 + t3547 * (t3365 * t3492 + t3366 * t3491 + t3367 * t3490); (t3551 + t3549 + t3496) * MDP(2) + (t3453 * t3496 + t3455 * t3551 + t3457 * t3549) * MDP(9) + (-t3447 * t3496 - t3449 * t3551 - t3451 * t3549) * MDP(10) + (t3425 * t3496 + t3426 * t3551 + t3427 * t3549) * MDP(11) + (-t3422 * t3496 - t3423 * t3551 - t3424 * t3549) * MDP(12) + ((t3458 * t3383 - (t3452 * t3443 + t3458 * t3540 + t3420) * t3544) * t3439 + (t3456 * t3385 - (t3450 * t3442 + t3456 * t3541 + t3419) * t3396) * t3438 + (t3454 * t3384 - (t3448 * t3441 + t3454 * t3542 + t3418) * t3395) * t3437) * MDP(14) - g(3) * MDP(15) + t3547 * (t3464 * t3508 + t3465 * t3509 + t3466 * t3510);];
taugX  = t1;
