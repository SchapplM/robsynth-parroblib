% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR6V1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR6V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:37
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR6V1G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:37:55
% EndTime: 2020-08-06 18:37:57
% DurationCPUTime: 2.13s
% Computational Cost: add. (912->194), mult. (1077->300), div. (111->10), fcn. (948->53), ass. (0->152)
t3450 = sin(pkin(7));
t3468 = -pkin(6) - pkin(5);
t3410 = t3450 * t3468 - pkin(1);
t3456 = sin(qJ(1,3));
t3401 = t3410 * t3456;
t3451 = cos(pkin(7));
t3462 = cos(qJ(1,3));
t3508 = t3462 * t3468;
t3520 = t3450 * t3462;
t3555 = t3401 - pkin(2) * t3520 - (pkin(2) * t3456 + t3508) * t3451;
t3458 = sin(qJ(1,2));
t3402 = t3410 * t3458;
t3464 = cos(qJ(1,2));
t3507 = t3464 * t3468;
t3519 = t3450 * t3464;
t3554 = t3402 - pkin(2) * t3519 - (pkin(2) * t3458 + t3507) * t3451;
t3460 = sin(qJ(1,1));
t3403 = t3410 * t3460;
t3466 = cos(qJ(1,1));
t3506 = t3466 * t3468;
t3518 = t3450 * t3466;
t3553 = t3403 - pkin(2) * t3518 - (pkin(2) * t3460 + t3506) * t3451;
t3539 = MDP(4) * pkin(1);
t3552 = MDP(2) / 0.2e1 + t3539 / 0.2e1;
t3550 = MDP(3) / 0.2e1;
t3549 = MDP(10) / 0.2e1;
t3548 = MDP(11) / 0.2e1;
t3543 = -0.2e1 * pkin(1);
t3542 = -0.2e1 * pkin(2);
t3541 = 0.2e1 * pkin(2);
t3540 = 0.2e1 * t3468;
t3461 = cos(qJ(3,3));
t3419 = t3461 * pkin(3) + pkin(2);
t3463 = cos(qJ(3,2));
t3420 = t3463 * pkin(3) + pkin(2);
t3465 = cos(qJ(3,1));
t3421 = t3465 * pkin(3) + pkin(2);
t3441 = qJ(1,3) + pkin(7);
t3422 = sin(t3441);
t3425 = cos(t3441);
t3455 = sin(qJ(3,3));
t3500 = pkin(7) + qJ(3,3);
t3503 = -pkin(7) + qJ(3,3);
t3538 = (t3422 * t3540 + t3462 * t3543 + t3425 * t3542 + (-cos(qJ(1,3) - t3503) - cos(qJ(1,3) + t3500)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,3)) + t3455 * t3541 + (sin(t3500) + sin(t3503)) * pkin(1));
t3442 = qJ(1,2) + pkin(7);
t3423 = sin(t3442);
t3426 = cos(t3442);
t3457 = sin(qJ(3,2));
t3501 = pkin(7) + qJ(3,2);
t3504 = -pkin(7) + qJ(3,2);
t3537 = (t3423 * t3540 + t3464 * t3543 + t3426 * t3542 + (-cos(qJ(1,2) - t3504) - cos(qJ(1,2) + t3501)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,2)) + t3457 * t3541 + (sin(t3501) + sin(t3504)) * pkin(1));
t3443 = qJ(1,1) + pkin(7);
t3424 = sin(t3443);
t3427 = cos(t3443);
t3459 = sin(qJ(3,1));
t3502 = pkin(7) + qJ(3,1);
t3505 = -pkin(7) + qJ(3,1);
t3536 = (t3424 * t3540 + t3466 * t3543 + t3427 * t3542 + (-cos(qJ(1,1) - t3505) - cos(qJ(1,1) + t3502)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,1)) + t3459 * t3541 + (sin(t3502) + sin(t3505)) * pkin(1));
t3452 = legFrame(3,2);
t3429 = sin(t3452);
t3432 = cos(t3452);
t3398 = g(1) * t3432 - g(2) * t3429;
t3428 = t3451 * pkin(1);
t3404 = 0.1e1 / (t3428 + t3419);
t3535 = (-g(3) * t3422 + t3398 * t3425) * t3404;
t3453 = legFrame(2,2);
t3430 = sin(t3453);
t3433 = cos(t3453);
t3399 = g(1) * t3433 - g(2) * t3430;
t3405 = 0.1e1 / (t3428 + t3420);
t3534 = (-g(3) * t3423 + t3399 * t3426) * t3405;
t3454 = legFrame(1,2);
t3431 = sin(t3454);
t3434 = cos(t3454);
t3400 = g(1) * t3434 - g(2) * t3431;
t3406 = 0.1e1 / (t3428 + t3421);
t3533 = (-g(3) * t3424 + t3400 * t3427) * t3406;
t3412 = t3452 + t3441;
t3413 = -t3452 + t3441;
t3386 = -sin(t3412) + sin(t3413);
t3532 = t3386 * t3404;
t3414 = t3453 + t3442;
t3415 = -t3453 + t3442;
t3387 = -sin(t3414) + sin(t3415);
t3531 = t3387 * t3405;
t3416 = t3454 + t3443;
t3417 = -t3454 + t3443;
t3388 = -sin(t3416) + sin(t3417);
t3530 = t3388 * t3406;
t3389 = cos(t3413) + cos(t3412);
t3529 = t3389 * t3404;
t3390 = cos(t3415) + cos(t3414);
t3528 = t3390 * t3405;
t3391 = cos(t3417) + cos(t3416);
t3527 = t3391 * t3406;
t3526 = t3404 * t3422;
t3525 = t3404 / t3455;
t3524 = t3405 * t3423;
t3523 = t3405 / t3457;
t3522 = t3406 * t3424;
t3521 = t3406 / t3459;
t3517 = t3455 * t3429;
t3516 = t3455 * t3432;
t3515 = t3456 * t3419;
t3514 = t3457 * t3430;
t3513 = t3457 * t3433;
t3512 = t3458 * t3420;
t3511 = t3459 * t3431;
t3510 = t3459 * t3434;
t3509 = t3460 * t3421;
t3499 = (t3456 * t3451 + t3520) * t3461 ^ 2 * pkin(3);
t3498 = (t3458 * t3451 + t3519) * t3463 ^ 2 * pkin(3);
t3497 = (t3460 * t3451 + t3518) * t3465 ^ 2 * pkin(3);
t3496 = ((t3508 + t3515) * t3451 - t3401 + t3419 * t3520) * t3525;
t3495 = ((t3507 + t3512) * t3451 - t3402 + t3420 * t3519) * t3523;
t3494 = ((t3506 + t3509) * t3451 - t3403 + t3421 * t3518) * t3521;
t3493 = t3455 * t3535;
t3492 = t3461 * t3535;
t3491 = t3457 * t3534;
t3490 = t3463 * t3534;
t3489 = t3459 * t3533;
t3488 = t3465 * t3533;
t3395 = g(1) * t3429 + g(2) * t3432;
t3487 = t3395 * t3525;
t3396 = g(1) * t3430 + g(2) * t3433;
t3486 = t3396 * t3523;
t3397 = g(1) * t3431 + g(2) * t3434;
t3485 = t3397 * t3521;
t3484 = g(3) * t3456 - t3398 * t3462;
t3483 = g(3) * t3458 - t3399 * t3464;
t3482 = g(3) * t3460 - t3400 * t3466;
t3481 = t3429 * t3496;
t3480 = t3432 * t3496;
t3479 = t3430 * t3495;
t3478 = t3433 * t3495;
t3477 = t3431 * t3494;
t3476 = t3434 * t3494;
t3475 = g(3) * t3425 + t3398 * t3422;
t3474 = g(3) * t3426 + t3399 * t3423;
t3473 = g(3) * t3427 + t3400 * t3424;
t3469 = 0.1e1 / pkin(3);
t3418 = t3428 + pkin(2);
t3385 = g(3) * t3466 + t3400 * t3460;
t3384 = g(3) * t3464 + t3399 * t3458;
t3383 = g(3) * t3462 + t3398 * t3456;
t3364 = t3397 * t3459 + t3473 * t3465;
t3363 = -t3397 * t3465 + t3473 * t3459;
t3362 = t3396 * t3457 + t3474 * t3463;
t3361 = -t3396 * t3463 + t3474 * t3457;
t3360 = t3395 * t3455 + t3475 * t3461;
t3359 = -t3395 * t3461 + t3475 * t3455;
t1 = [(-(t3434 * t3497 + (pkin(3) * t3511 - t3434 * t3553) * t3465 + t3418 * t3511) * t3485 - (t3433 * t3498 + (pkin(3) * t3514 - t3433 * t3554) * t3463 + t3418 * t3514) * t3486 - (t3432 * t3499 + (pkin(3) * t3517 - t3432 * t3555) * t3461 + t3418 * t3517) * t3487) * MDP(4) - g(1) * MDP(12) + ((-t3359 * t3480 - t3361 * t3478 - t3363 * t3476) * MDP(10) + (-t3360 * t3480 - t3362 * t3478 - t3364 * t3476) * MDP(11)) * t3469 + (t3383 * t3529 + t3384 * t3528 + t3385 * t3527) * t3550 + (-t3389 * t3492 - t3390 * t3490 - t3391 * t3488) * t3549 + (t3389 * t3493 + t3390 * t3491 + t3391 * t3489) * t3548 + t3552 * (t3482 * t3527 + t3483 * t3528 + t3484 * t3529); (-(-t3431 * t3497 + (pkin(3) * t3510 + t3431 * t3553) * t3465 + t3418 * t3510) * t3485 - (-t3430 * t3498 + (pkin(3) * t3513 + t3430 * t3554) * t3463 + t3418 * t3513) * t3486 - (-t3429 * t3499 + (pkin(3) * t3516 + t3429 * t3555) * t3461 + t3418 * t3516) * t3487) * MDP(4) - g(2) * MDP(12) + ((t3359 * t3481 + t3361 * t3479 + t3363 * t3477) * MDP(10) + (t3360 * t3481 + t3362 * t3479 + t3364 * t3477) * MDP(11)) * t3469 + (t3383 * t3532 + t3384 * t3531 + t3385 * t3530) * t3550 + (-t3386 * t3492 - t3387 * t3490 - t3388 * t3488) * t3549 + (t3386 * t3493 + t3387 * t3491 + t3388 * t3489) * t3548 + t3552 * (t3482 * t3530 + t3483 * t3531 + t3484 * t3532); (-t3383 * t3526 - t3384 * t3524 - t3385 * t3522) * MDP(3) + (-((t3421 * t3466 - t3460 * t3468) * t3451 - t3410 * t3466 - t3450 * t3509) * t3465 * t3485 - ((t3420 * t3464 - t3458 * t3468) * t3451 - t3410 * t3464 - t3450 * t3512) * t3463 * t3486 - ((t3419 * t3462 - t3456 * t3468) * t3451 - t3410 * t3462 - t3450 * t3515) * t3461 * t3487) * MDP(4) + (t3422 * t3492 + t3423 * t3490 + t3424 * t3488) * MDP(10) + (-t3422 * t3493 - t3423 * t3491 - t3424 * t3489) * MDP(11) - g(3) * MDP(12) + ((t3359 * t3538 + t3361 * t3537 + t3363 * t3536) * MDP(10) + (t3360 * t3538 + t3362 * t3537 + t3364 * t3536) * MDP(11)) * t3469 + (MDP(2) + t3539) * (-t3482 * t3522 - t3483 * t3524 - t3484 * t3526);];
taugX  = t1;
