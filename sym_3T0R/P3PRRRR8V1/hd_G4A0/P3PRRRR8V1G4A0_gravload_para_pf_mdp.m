% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR8V1G4A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR8V1G4A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:27
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V1G4A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:27:15
% EndTime: 2020-08-06 17:27:19
% DurationCPUTime: 3.68s
% Computational Cost: add. (1869->264), mult. (4590->549), div. (126->7), fcn. (5253->34), ass. (0->237)
t3648 = legFrame(3,2);
t3629 = sin(t3648);
t3632 = cos(t3648);
t3645 = legFrame(3,1);
t3620 = sin(t3645);
t3626 = cos(t3645);
t3772 = -g(2) * t3620 + g(3) * t3626;
t3557 = t3629 * g(1) + t3772 * t3632;
t3649 = legFrame(2,2);
t3630 = sin(t3649);
t3633 = cos(t3649);
t3646 = legFrame(2,1);
t3621 = sin(t3646);
t3627 = cos(t3646);
t3773 = -g(2) * t3621 + g(3) * t3627;
t3558 = t3630 * g(1) + t3773 * t3633;
t3650 = legFrame(1,2);
t3631 = sin(t3650);
t3634 = cos(t3650);
t3647 = legFrame(1,1);
t3622 = sin(t3647);
t3628 = cos(t3647);
t3774 = -g(2) * t3622 + g(3) * t3628;
t3559 = t3631 * g(1) + t3774 * t3634;
t3661 = cos(qJ(3,1));
t3637 = 0.1e1 / t3661;
t3662 = cos(qJ(2,1));
t3656 = sin(qJ(2,1));
t3706 = t3656 * t3661;
t3601 = pkin(2) * t3706 - t3662 * pkin(5);
t3639 = sin(pkin(3));
t3641 = cos(pkin(3));
t3655 = sin(qJ(3,1));
t3722 = t3641 * t3655;
t3571 = pkin(2) * t3722 + t3601 * t3639;
t3763 = 0.1e1 / t3571;
t3777 = t3763 * t3637;
t3659 = cos(qJ(3,2));
t3636 = 0.1e1 / t3659;
t3660 = cos(qJ(2,2));
t3654 = sin(qJ(2,2));
t3710 = t3654 * t3659;
t3600 = pkin(2) * t3710 - t3660 * pkin(5);
t3653 = sin(qJ(3,2));
t3724 = t3641 * t3653;
t3570 = pkin(2) * t3724 + t3600 * t3639;
t3764 = 0.1e1 / t3570;
t3776 = t3764 * t3636;
t3657 = cos(qJ(3,3));
t3635 = 0.1e1 / t3657;
t3658 = cos(qJ(2,3));
t3652 = sin(qJ(2,3));
t3714 = t3652 * t3657;
t3599 = pkin(2) * t3714 - t3658 * pkin(5);
t3651 = sin(qJ(3,3));
t3726 = t3641 * t3651;
t3569 = pkin(2) * t3726 + t3599 * t3639;
t3765 = 0.1e1 / t3569;
t3775 = t3765 * t3635;
t3593 = t3626 * g(2) + t3620 * g(3);
t3642 = legFrame(3,3);
t3617 = sin(t3642);
t3623 = cos(t3642);
t3684 = t3632 * g(1) - t3629 * t3772;
t3771 = t3617 * t3593 + t3623 * t3684;
t3594 = t3627 * g(2) + t3621 * g(3);
t3643 = legFrame(2,3);
t3618 = sin(t3643);
t3624 = cos(t3643);
t3683 = t3633 * g(1) - t3630 * t3773;
t3770 = t3618 * t3594 + t3624 * t3683;
t3595 = t3628 * g(2) + t3622 * g(3);
t3644 = legFrame(1,3);
t3619 = sin(t3644);
t3625 = cos(t3644);
t3682 = t3634 * g(1) - t3631 * t3774;
t3769 = t3619 * t3595 + t3625 * t3682;
t3768 = -t3623 * t3593 + t3617 * t3684;
t3767 = -t3624 * t3594 + t3618 * t3683;
t3766 = -t3625 * t3595 + t3619 * t3682;
t3762 = pkin(2) * t3657;
t3761 = pkin(2) * t3659;
t3760 = pkin(2) * t3661;
t3638 = sin(pkin(6));
t3640 = cos(pkin(6));
t3533 = t3638 * t3593 + t3640 * t3684;
t3536 = -t3593 * t3640 + t3638 * t3684;
t3729 = t3639 * t3557;
t3476 = (-t3729 + (t3533 * t3617 + t3536 * t3623) * t3641) * t3658 + t3652 * (t3533 * t3623 - t3536 * t3617);
t3750 = t3476 * t3765;
t3534 = t3638 * t3594 + t3640 * t3683;
t3537 = -t3594 * t3640 + t3638 * t3683;
t3728 = t3639 * t3558;
t3477 = (-t3728 + (t3534 * t3618 + t3537 * t3624) * t3641) * t3660 + t3654 * (t3534 * t3624 - t3537 * t3618);
t3749 = t3477 * t3764;
t3535 = t3638 * t3595 + t3640 * t3682;
t3538 = -t3595 * t3640 + t3638 * t3682;
t3727 = t3639 * t3559;
t3478 = (-t3727 + (t3535 * t3619 + t3538 * t3625) * t3641) * t3662 + t3656 * (t3535 * t3625 - t3538 * t3619);
t3748 = t3478 * t3763;
t3747 = t3557 * t3765;
t3746 = t3558 * t3764;
t3745 = t3559 * t3763;
t3738 = t3620 * t3629;
t3737 = t3621 * t3630;
t3736 = t3622 * t3631;
t3732 = t3626 * t3629;
t3731 = t3627 * t3630;
t3730 = t3628 * t3631;
t3725 = t3641 * t3652;
t3723 = t3641 * t3654;
t3721 = t3641 * t3656;
t3720 = t3641 * t3658;
t3719 = t3641 * t3660;
t3718 = t3641 * t3662;
t3717 = t3651 * t3639;
t3716 = t3651 * t3652;
t3715 = t3651 * t3658;
t3713 = t3653 * t3639;
t3712 = t3653 * t3654;
t3711 = t3653 * t3660;
t3709 = t3655 * t3639;
t3708 = t3655 * t3656;
t3707 = t3655 * t3662;
t3705 = t3657 * t3639;
t3704 = t3657 * t3658;
t3703 = t3659 * t3639;
t3702 = t3659 * t3660;
t3701 = t3661 * t3639;
t3700 = t3661 * t3662;
t3587 = t3641 * t3716 + t3705;
t3485 = t3771 * (-t3587 * t3638 + t3640 * t3715) - t3768 * (t3587 * t3640 + t3638 * t3715) - t3557 * (-t3639 * t3716 + t3641 * t3657);
t3699 = t3485 * t3775;
t3588 = t3641 * t3712 + t3703;
t3486 = t3770 * (-t3588 * t3638 + t3640 * t3711) - t3767 * (t3588 * t3640 + t3638 * t3711) - t3558 * (-t3639 * t3712 + t3641 * t3659);
t3698 = t3486 * t3776;
t3589 = t3641 * t3708 + t3701;
t3487 = t3769 * (-t3589 * t3638 + t3640 * t3707) - (t3589 * t3640 + t3638 * t3707) * t3766 - t3559 * (-t3639 * t3708 + t3641 * t3661);
t3697 = t3487 * t3777;
t3590 = t3641 * t3714 - t3717;
t3488 = (-t3590 * t3638 + t3640 * t3704) * t3771 - (t3590 * t3640 + t3638 * t3704) * t3768 + t3557 * (t3652 * t3705 + t3726);
t3696 = t3488 * t3775;
t3591 = t3641 * t3710 - t3713;
t3489 = (-t3591 * t3638 + t3640 * t3702) * t3770 - (t3591 * t3640 + t3638 * t3702) * t3767 + t3558 * (t3654 * t3703 + t3724);
t3695 = t3489 * t3776;
t3592 = t3641 * t3706 - t3709;
t3490 = t3769 * (-t3592 * t3638 + t3640 * t3700) - (t3592 * t3640 + t3638 * t3700) * t3766 + (t3656 * t3701 + t3722) * t3559;
t3694 = t3490 * t3777;
t3575 = t3617 * t3640 + t3623 * t3638;
t3578 = -t3617 * t3638 + t3623 * t3640;
t3522 = -t3620 * t3575 + t3578 * t3732;
t3582 = t3638 * t3658 + t3640 * t3725;
t3585 = -t3638 * t3725 + t3640 * t3658;
t3551 = t3582 * t3623 + t3617 * t3585;
t3668 = t3617 * t3582 - t3585 * t3623;
t3497 = (t3551 * t3732 - t3668 * t3620) * t3651 + t3522 * t3705;
t3693 = t3497 * t3775;
t3576 = t3618 * t3640 + t3624 * t3638;
t3579 = -t3618 * t3638 + t3624 * t3640;
t3524 = -t3621 * t3576 + t3579 * t3731;
t3583 = t3638 * t3660 + t3640 * t3723;
t3586 = -t3638 * t3723 + t3640 * t3660;
t3552 = t3583 * t3624 + t3618 * t3586;
t3667 = t3618 * t3583 - t3586 * t3624;
t3498 = (t3552 * t3731 - t3667 * t3621) * t3653 + t3524 * t3703;
t3692 = t3498 * t3776;
t3577 = t3619 * t3640 + t3625 * t3638;
t3580 = -t3619 * t3638 + t3625 * t3640;
t3526 = -t3622 * t3577 + t3580 * t3730;
t3581 = t3638 * t3721 - t3640 * t3662;
t3584 = t3638 * t3662 + t3640 * t3721;
t3553 = -t3619 * t3581 + t3584 * t3625;
t3669 = t3581 * t3625 + t3619 * t3584;
t3499 = (t3553 * t3730 - t3669 * t3622) * t3655 + t3526 * t3701;
t3691 = t3499 * t3777;
t3527 = t3575 * t3626 + t3578 * t3738;
t3500 = (-t3551 * t3738 - t3668 * t3626) * t3651 - t3527 * t3705;
t3690 = t3500 * t3775;
t3528 = t3576 * t3627 + t3579 * t3737;
t3501 = (-t3552 * t3737 - t3667 * t3627) * t3653 - t3528 * t3703;
t3689 = t3501 * t3776;
t3529 = t3577 * t3628 + t3580 * t3736;
t3502 = (-t3553 * t3736 - t3669 * t3628) * t3655 - t3529 * t3701;
t3688 = t3502 * t3777;
t3687 = (t3551 * t3651 + t3578 * t3705) * t3765 * t3632;
t3686 = (t3552 * t3653 + t3579 * t3703) * t3764 * t3633;
t3685 = (t3553 * t3655 + t3580 * t3701) * t3763 * t3634;
t3681 = t3476 * t3651 * t3775;
t3680 = t3477 * t3653 * t3776;
t3679 = t3478 * t3655 * t3777;
t3678 = ((-t3652 * t3575 + t3578 * t3720) * t3762 + pkin(5) * (t3658 * t3575 + t3578 * t3725)) * t3632 * t3775;
t3677 = ((-t3654 * t3576 + t3579 * t3719) * t3761 + pkin(5) * (t3660 * t3576 + t3579 * t3723)) * t3633 * t3776;
t3676 = ((-t3656 * t3577 + t3580 * t3718) * t3760 + pkin(5) * (t3662 * t3577 + t3580 * t3721)) * t3634 * t3777;
t3675 = t3635 * t3687;
t3674 = t3476 * t3687;
t3673 = t3636 * t3686;
t3672 = t3477 * t3686;
t3671 = t3637 * t3685;
t3670 = t3478 * t3685;
t3666 = pkin(2) * t3717 - t3599 * t3641;
t3665 = pkin(2) * t3713 - t3600 * t3641;
t3664 = pkin(2) * t3709 - t3601 * t3641;
t3663 = 0.1e1 / pkin(2);
t3604 = pkin(2) * t3700 + pkin(5) * t3656;
t3603 = pkin(2) * t3702 + pkin(5) * t3654;
t3602 = pkin(2) * t3704 + pkin(5) * t3652;
t3550 = -t3638 * t3604 + t3664 * t3640;
t3549 = -t3638 * t3603 + t3665 * t3640;
t3548 = -t3638 * t3602 + t3666 * t3640;
t3547 = t3604 * t3640 + t3664 * t3638;
t3546 = t3603 * t3640 + t3665 * t3638;
t3545 = t3602 * t3640 + t3666 * t3638;
t3532 = t3577 * t3736 - t3580 * t3628;
t3531 = t3576 * t3737 - t3579 * t3627;
t3530 = t3575 * t3738 - t3578 * t3626;
t3525 = t3577 * t3730 + t3622 * t3580;
t3523 = t3576 * t3731 + t3621 * t3579;
t3521 = t3575 * t3732 + t3620 * t3578;
t3508 = -t3547 * t3619 + t3550 * t3625;
t3507 = -t3546 * t3618 + t3549 * t3624;
t3506 = -t3545 * t3617 + t3548 * t3623;
t3505 = -t3634 * t3571 + (t3547 * t3625 + t3550 * t3619) * t3631;
t3504 = -t3633 * t3570 + (t3546 * t3624 + t3549 * t3618) * t3630;
t3503 = -t3632 * t3569 + (t3545 * t3623 + t3548 * t3617) * t3629;
t3496 = t3769 * (t3638 * t3718 + t3640 * t3656) + t3766 * (-t3638 * t3656 + t3640 * t3718) - t3662 * t3727;
t3495 = t3770 * (t3638 * t3719 + t3640 * t3654) + t3767 * (-t3638 * t3654 + t3640 * t3719) - t3660 * t3728;
t3494 = -t3581 * t3769 - t3584 * t3766 + t3656 * t3727;
t3493 = -t3583 * t3767 + t3586 * t3770 + t3654 * t3728;
t3492 = -t3582 * t3768 + t3585 * t3771 + t3652 * t3729;
t3491 = (t3638 * t3720 + t3640 * t3652) * t3771 + t3768 * (-t3638 * t3652 + t3640 * t3720) - t3658 * t3729;
t3484 = -(t3529 * t3718 - t3532 * t3656) * t3760 - pkin(5) * (t3529 * t3721 + t3662 * t3532);
t3483 = (-t3656 * t3525 + t3526 * t3718) * t3760 + pkin(5) * (t3662 * t3525 + t3526 * t3721);
t3482 = -(t3528 * t3719 - t3531 * t3654) * t3761 - pkin(5) * (t3528 * t3723 + t3660 * t3531);
t3481 = (-t3654 * t3523 + t3524 * t3719) * t3761 + pkin(5) * (t3660 * t3523 + t3524 * t3723);
t3480 = -(t3527 * t3720 - t3530 * t3652) * t3762 - pkin(5) * (t3527 * t3725 + t3658 * t3530);
t3479 = (-t3652 * t3521 + t3522 * t3720) * t3762 + pkin(5) * (t3658 * t3521 + t3522 * t3725);
t1 = [(-((t3664 * t3577 + t3580 * t3604) * t3634 + t3631 * t3571) * t3745 - ((t3665 * t3576 + t3579 * t3603) * t3633 + t3630 * t3570) * t3746 - ((t3666 * t3575 + t3578 * t3602) * t3632 + t3629 * t3569) * t3747) * MDP(1) + (-t3491 * t3675 - t3495 * t3673 - t3496 * t3671) * MDP(3) + (-t3492 * t3675 - t3493 * t3673 - t3494 * t3671) * MDP(4) + (-t3670 - t3672 - t3674) * MDP(10) + (t3635 * t3651 * t3674 + t3636 * t3653 * t3672 + t3637 * t3655 * t3670) * MDP(11) - g(1) * MDP(12) + ((-t3485 * t3678 - t3486 * t3677 - t3487 * t3676) * MDP(10) + (-t3488 * t3678 - t3489 * t3677 - t3490 * t3676) * MDP(11)) * t3663; (-(t3505 * t3622 - t3628 * t3508) * t3745 - (t3504 * t3621 - t3627 * t3507) * t3746 - (t3503 * t3620 - t3626 * t3506) * t3747) * MDP(1) + (t3491 * t3690 + t3495 * t3689 + t3496 * t3688) * MDP(3) + (t3492 * t3690 + t3493 * t3689 + t3494 * t3688) * MDP(4) + (t3500 * t3750 + t3501 * t3749 + t3502 * t3748) * MDP(10) + (-t3500 * t3681 - t3501 * t3680 - t3502 * t3679) * MDP(11) - g(2) * MDP(12) + ((t3480 * t3699 + t3482 * t3698 + t3484 * t3697) * MDP(10) + (t3480 * t3696 + t3482 * t3695 + t3484 * t3694) * MDP(11)) * t3663; (-(-t3505 * t3628 - t3622 * t3508) * t3745 - (-t3504 * t3627 - t3621 * t3507) * t3746 - (-t3503 * t3626 - t3620 * t3506) * t3747) * MDP(1) + (t3491 * t3693 + t3495 * t3692 + t3496 * t3691) * MDP(3) + (t3492 * t3693 + t3493 * t3692 + t3494 * t3691) * MDP(4) + (t3497 * t3750 + t3498 * t3749 + t3499 * t3748) * MDP(10) + (-t3497 * t3681 - t3498 * t3680 - t3499 * t3679) * MDP(11) - g(3) * MDP(12) + ((t3479 * t3699 + t3481 * t3698 + t3483 * t3697) * MDP(10) + (t3479 * t3696 + t3481 * t3695 + t3483 * t3694) * MDP(11)) * t3663;];
taugX  = t1;
