% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR8V2G2A0
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [13x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR8V2G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 21:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V2G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(13,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_mdp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:13:36
% EndTime: 2020-08-06 21:13:38
% DurationCPUTime: 1.60s
% Computational Cost: add. (927->191), mult. (1608->328), div. (123->9), fcn. (1374->23), ass. (0->168)
t3646 = MDP(3) - MDP(11);
t3544 = sin(qJ(1,3));
t3550 = cos(qJ(1,3));
t3540 = legFrame(3,2);
t3515 = sin(t3540);
t3518 = cos(t3540);
t3568 = t3518 * g(1) - t3515 * g(2);
t3479 = g(3) * t3544 - t3568 * t3550;
t3546 = sin(qJ(1,2));
t3552 = cos(qJ(1,2));
t3541 = legFrame(2,2);
t3516 = sin(t3541);
t3519 = cos(t3541);
t3567 = t3519 * g(1) - t3516 * g(2);
t3480 = g(3) * t3546 - t3567 * t3552;
t3548 = sin(qJ(1,1));
t3554 = cos(qJ(1,1));
t3542 = legFrame(1,2);
t3517 = sin(t3542);
t3520 = cos(t3542);
t3566 = t3520 * g(1) - t3517 * g(2);
t3481 = g(3) * t3548 - t3566 * t3554;
t3549 = cos(qJ(2,3));
t3645 = 0.2e1 * t3549 ^ 2;
t3551 = cos(qJ(2,2));
t3644 = 0.2e1 * t3551 ^ 2;
t3553 = cos(qJ(2,1));
t3643 = 0.2e1 * t3553 ^ 2;
t3539 = -qJ(3,1) - pkin(5);
t3538 = -qJ(3,2) - pkin(5);
t3537 = -qJ(3,3) - pkin(5);
t3642 = MDP(12) * pkin(2);
t3528 = pkin(6) - t3537;
t3577 = pkin(1) * t3544 - t3550 * t3528;
t3536 = cos(pkin(7));
t3531 = t3536 ^ 2;
t3602 = pkin(3) * (t3531 - 0.1e1);
t3535 = sin(pkin(7));
t3543 = sin(qJ(2,3));
t3606 = t3535 * t3543;
t3641 = pkin(3) * (t3544 * t3602 + t3577 * t3606);
t3529 = pkin(6) - t3538;
t3576 = pkin(1) * t3546 - t3552 * t3529;
t3545 = sin(qJ(2,2));
t3605 = t3535 * t3545;
t3640 = pkin(3) * (t3546 * t3602 + t3576 * t3605);
t3530 = pkin(6) - t3539;
t3575 = pkin(1) * t3548 - t3554 * t3530;
t3547 = sin(qJ(2,1));
t3604 = t3535 * t3547;
t3639 = pkin(3) * (t3548 * t3602 + t3575 * t3604);
t3638 = pkin(3) * cos(qJ(2,3) + pkin(7));
t3637 = pkin(3) * cos(qJ(2,2) + pkin(7));
t3636 = pkin(3) * cos(qJ(2,1) + pkin(7));
t3632 = t3535 * pkin(3);
t3631 = t3536 * pkin(3);
t3504 = pkin(2) + t3631;
t3595 = pkin(3) * t3606;
t3485 = 0.1e1 / (t3504 * t3549 - t3595);
t3512 = 0.1e1 / t3528;
t3630 = t3485 * t3512;
t3594 = pkin(3) * t3605;
t3486 = 0.1e1 / (t3504 * t3551 - t3594);
t3513 = 0.1e1 / t3529;
t3629 = t3486 * t3513;
t3593 = pkin(3) * t3604;
t3487 = 0.1e1 / (t3504 * t3553 - t3593);
t3514 = 0.1e1 / t3530;
t3628 = t3487 * t3514;
t3521 = t3549 * pkin(2);
t3495 = 0.1e1 / (t3521 + t3638);
t3627 = t3495 * t3515;
t3626 = t3495 * t3518;
t3522 = t3551 * pkin(2);
t3496 = 0.1e1 / (t3522 + t3637);
t3625 = t3496 * t3516;
t3624 = t3496 * t3519;
t3523 = t3553 * pkin(2);
t3497 = 0.1e1 / (t3523 + t3636);
t3623 = t3497 * t3517;
t3622 = t3497 * t3520;
t3621 = t3504 * t3518;
t3620 = t3504 * t3519;
t3619 = t3504 * t3520;
t3618 = t3512 * t3550;
t3617 = t3513 * t3552;
t3616 = t3514 * t3554;
t3615 = t3515 * t3504;
t3614 = t3515 * t3544;
t3613 = t3516 * t3504;
t3612 = t3516 * t3546;
t3611 = t3517 * t3504;
t3610 = t3517 * t3548;
t3609 = t3518 * t3544;
t3608 = t3519 * t3546;
t3607 = t3520 * t3548;
t3603 = pkin(2) * t3631;
t3601 = t3518 * t3632;
t3600 = t3519 * t3632;
t3599 = t3520 * t3632;
t3598 = t3515 * t3632;
t3597 = t3516 * t3632;
t3596 = t3517 * t3632;
t3592 = t3479 * t3512 * t3549;
t3591 = t3480 * t3513 * t3551;
t3590 = t3481 * t3514 * t3553;
t3589 = t3479 * t3630;
t3588 = t3480 * t3629;
t3587 = t3481 * t3628;
t3564 = g(3) * t3550 + t3568 * t3544;
t3586 = t3564 * t3630;
t3563 = g(3) * t3552 + t3567 * t3546;
t3585 = t3563 * t3629;
t3562 = g(3) * t3554 + t3566 * t3548;
t3584 = t3562 * t3628;
t3505 = t3521 + pkin(1);
t3502 = t3505 * t3550;
t3457 = g(3) * (t3544 * t3505 + t3537 * t3550) - t3568 * (-t3544 * t3537 + t3502);
t3583 = t3457 * t3630;
t3506 = t3522 + pkin(1);
t3503 = t3506 * t3552;
t3455 = (t3546 * t3506 + t3538 * t3552) * g(3) - t3567 * (-t3546 * t3538 + t3503);
t3582 = t3455 * t3629;
t3507 = t3523 + pkin(1);
t3501 = t3554 * t3507;
t3456 = (t3548 * t3507 + t3539 * t3554) * g(3) - t3566 * (-t3548 * t3539 + t3501);
t3581 = t3456 * t3628;
t3580 = t3479 * t3618;
t3579 = t3480 * t3617;
t3578 = t3481 * t3616;
t3574 = t3485 * t3592;
t3573 = t3486 * t3591;
t3572 = t3487 * t3590;
t3571 = t3543 * t3589;
t3570 = t3545 * t3588;
t3569 = t3547 * t3587;
t3555 = pkin(3) ^ 2;
t3556 = pkin(2) ^ 2;
t3565 = 0.2e1 * t3531 * t3555 - t3555 + t3556 + 0.2e1 * t3603;
t3492 = t3515 * g(1) + t3518 * g(2);
t3458 = -t3492 * t3549 + t3543 * t3564;
t3493 = t3516 * g(1) + t3519 * g(2);
t3460 = -t3493 * t3551 + t3545 * t3563;
t3494 = t3517 * g(1) + t3520 * g(2);
t3462 = -t3494 * t3553 + t3547 * t3562;
t3561 = t3458 * t3627 + t3460 * t3625 + t3462 * t3623;
t3560 = t3458 * t3626 + t3460 * t3624 + t3462 * t3622;
t3508 = pkin(1) * t3632;
t3500 = pkin(1) * t3547 - t3632;
t3499 = pkin(1) * t3545 - t3632;
t3498 = pkin(1) * t3543 - t3632;
t3491 = t3603 + t3556 / 0.2e1 + (t3531 - 0.1e1 / 0.2e1) * t3555;
t3472 = -0.2e1 * t3548 * t3593 + t3575;
t3471 = -0.2e1 * t3546 * t3594 + t3576;
t3470 = -0.2e1 * t3544 * t3595 + t3577;
t3469 = t3565 * t3547 + t3508;
t3468 = t3565 * t3545 + t3508;
t3467 = t3565 * t3543 + t3508;
t3463 = t3494 * t3547 + t3553 * t3562;
t3461 = t3493 * t3545 + t3551 * t3563;
t3459 = t3492 * t3543 + t3549 * t3564;
t3454 = (t3504 * t3607 + t3596) * t3553 + t3547 * (-t3548 * t3599 + t3611);
t3453 = (t3504 * t3608 + t3597) * t3551 + t3545 * (-t3546 * t3600 + t3613);
t3452 = (t3504 * t3609 + t3598) * t3549 + t3543 * (-t3544 * t3601 + t3615);
t3451 = (-t3504 * t3610 + t3599) * t3553 + (t3548 * t3596 + t3619) * t3547;
t3450 = (-t3504 * t3612 + t3600) * t3551 + (t3546 * t3597 + t3620) * t3545;
t3449 = (-t3504 * t3614 + t3601) * t3549 + (t3544 * t3598 + t3621) * t3543;
t1 = [(t3452 * t3589 + t3453 * t3588 + t3454 * t3587) * MDP(2) + (t3452 * t3574 + t3453 * t3573 + t3454 * t3572 + t3561) * MDP(9) + (-t3452 * t3571 - t3453 * t3570 - t3454 * t3569 + t3459 * t3627 + t3461 * t3625 + t3463 * t3623) * MDP(10) + (t3454 * t3581 - ((t3491 * t3607 + t3504 * t3596) * t3643 + (t3517 * t3469 + t3472 * t3619) * t3553 - t3520 * t3639 + t3500 * t3611) * t3587 + t3453 * t3582 - ((t3491 * t3608 + t3504 * t3597) * t3644 + (t3516 * t3468 + t3471 * t3620) * t3551 - t3519 * t3640 + t3499 * t3613) * t3588 + t3452 * t3583 - ((t3491 * t3609 + t3504 * t3598) * t3645 + (t3515 * t3467 + t3470 * t3621) * t3549 - t3518 * t3641 + t3498 * t3615) * t3589) * MDP(12) - g(1) * MDP(13) + t3561 * t3642 + t3646 * (t3452 * t3586 + t3453 * t3585 + t3454 * t3584); (t3449 * t3589 + t3450 * t3588 + t3451 * t3587) * MDP(2) + (t3449 * t3574 + t3450 * t3573 + t3451 * t3572 + t3560) * MDP(9) + (-t3449 * t3571 - t3450 * t3570 - t3451 * t3569 + t3459 * t3626 + t3461 * t3624 + t3463 * t3622) * MDP(10) + (t3451 * t3581 - ((-t3491 * t3610 + t3504 * t3599) * t3643 + (t3469 * t3520 - t3472 * t3611) * t3553 + t3517 * t3639 + t3500 * t3619) * t3587 + t3450 * t3582 - ((-t3491 * t3612 + t3504 * t3600) * t3644 + (t3468 * t3519 - t3471 * t3613) * t3551 + t3516 * t3640 + t3499 * t3620) * t3588 + t3449 * t3583 - ((-t3491 * t3614 + t3504 * t3601) * t3645 + (t3467 * t3518 - t3470 * t3615) * t3549 + t3515 * t3641 + t3498 * t3621) * t3589) * MDP(12) - g(2) * MDP(13) + t3560 * t3642 + t3646 * (t3449 * t3586 + t3450 * t3585 + t3451 * t3584); (t3578 + t3579 + t3580) * MDP(2) + (t3550 * t3592 + t3552 * t3591 + t3554 * t3590) * MDP(9) + (-t3543 * t3580 - t3545 * t3579 - t3547 * t3578) * MDP(10) + ((t3554 * t3456 - (t3548 * t3530 + t3554 * t3636 + t3501) * t3481) * t3514 + (t3552 * t3455 - (t3546 * t3529 + t3552 * t3637 + t3503) * t3480) * t3513 + (t3550 * t3457 - (t3544 * t3528 + t3550 * t3638 + t3502) * t3479) * t3512) * MDP(12) - g(3) * MDP(13) + t3646 * (t3562 * t3616 + t3563 * t3617 + t3564 * t3618);];
taugX  = t1;
