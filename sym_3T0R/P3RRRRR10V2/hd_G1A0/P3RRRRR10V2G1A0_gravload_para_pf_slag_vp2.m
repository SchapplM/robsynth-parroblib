% Calculate Gravitation load for parallel robot
% P3RRRRR10V2G1A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 00:48
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR10V2G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V2G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 23:49:29
% EndTime: 2020-08-06 23:49:34
% DurationCPUTime: 6.06s
% Computational Cost: add. (2301->448), mult. (3771->734), div. (36->16), fcn. (2877->26), ass. (0->270)
t3841 = m(2) + m(3);
t3680 = legFrame(1,3);
t3648 = sin(t3680);
t3651 = cos(t3680);
t3590 = -t3648 * g(1) + t3651 * g(2);
t3593 = t3651 * g(1) + t3648 * g(2);
t3689 = sin(qJ(1,1));
t3698 = cos(qJ(1,1));
t3541 = t3689 * t3590 + t3698 * t3593;
t3688 = sin(qJ(2,1));
t3697 = cos(qJ(2,1));
t3677 = cos(pkin(4));
t3725 = t3590 * t3698 - t3593 * t3689;
t3676 = sin(pkin(4));
t3835 = g(3) * t3676;
t3728 = t3725 * t3677 + t3835;
t3707 = t3697 * t3541 + t3728 * t3688;
t3679 = legFrame(2,3);
t3647 = sin(t3679);
t3650 = cos(t3679);
t3589 = -t3647 * g(1) + t3650 * g(2);
t3592 = t3650 * g(1) + t3647 * g(2);
t3686 = sin(qJ(1,2));
t3695 = cos(qJ(1,2));
t3540 = t3686 * t3589 + t3695 * t3592;
t3685 = sin(qJ(2,2));
t3694 = cos(qJ(2,2));
t3726 = t3589 * t3695 - t3592 * t3686;
t3729 = t3726 * t3677 + t3835;
t3708 = t3694 * t3540 + t3729 * t3685;
t3678 = legFrame(3,3);
t3646 = sin(t3678);
t3649 = cos(t3678);
t3588 = -t3646 * g(1) + t3649 * g(2);
t3591 = t3649 * g(1) + t3646 * g(2);
t3683 = sin(qJ(1,3));
t3692 = cos(qJ(1,3));
t3539 = t3683 * t3588 + t3692 * t3591;
t3682 = sin(qJ(2,3));
t3691 = cos(qJ(2,3));
t3727 = t3588 * t3692 - t3591 * t3683;
t3730 = t3727 * t3677 + t3835;
t3709 = t3691 * t3539 + t3730 * t3682;
t3834 = g(3) * t3677;
t3690 = cos(qJ(3,3));
t3660 = t3690 * pkin(3);
t3626 = t3660 + pkin(2);
t3693 = cos(qJ(3,2));
t3662 = t3693 * pkin(3);
t3628 = t3662 + pkin(2);
t3696 = cos(qJ(3,1));
t3664 = t3696 * pkin(3);
t3630 = t3664 + pkin(2);
t3700 = pkin(8) + pkin(7);
t3837 = pkin(1) * t3677;
t3836 = pkin(2) * t3677;
t3670 = t3690 ^ 2;
t3833 = t3670 * pkin(3);
t3672 = t3693 ^ 2;
t3832 = t3672 * pkin(3);
t3674 = t3696 ^ 2;
t3831 = t3674 * pkin(3);
t3681 = sin(qJ(3,3));
t3830 = t3681 * pkin(2);
t3657 = t3681 * pkin(3);
t3684 = sin(qJ(3,2));
t3829 = t3684 * pkin(2);
t3658 = t3684 * pkin(3);
t3687 = sin(qJ(3,1));
t3828 = t3687 * pkin(2);
t3659 = t3687 * pkin(3);
t3661 = t3691 * pkin(2);
t3663 = t3694 * pkin(2);
t3665 = t3697 * pkin(2);
t3667 = t3677 ^ 2;
t3827 = (t3667 - 0.1e1) * pkin(6);
t3826 = mrSges(3,2) * t3676;
t3825 = -0.2e1 * pkin(2) * pkin(3);
t3824 = -m(3) * pkin(2) - mrSges(2,1);
t3669 = pkin(2) - t3700;
t3668 = pkin(2) + t3700;
t3635 = t3700 * t3691;
t3600 = -t3682 * pkin(2) + t3635;
t3575 = t3600 * t3837;
t3621 = t3657 - pkin(6);
t3585 = -mrSges(3,1) * t3690 + mrSges(3,2) * t3681 + t3824;
t3616 = t3841 * pkin(1) + mrSges(1,1);
t3638 = pkin(7) * m(3) - mrSges(2,2) + mrSges(3,3);
t3718 = t3585 * t3691 - t3638 * t3682 - t3616;
t3762 = t3841 * pkin(6) + mrSges(2,3);
t3721 = (-t3682 * t3585 - t3638 * t3691) * t3677 - (mrSges(3,1) * t3681 + mrSges(3,2) * t3690 + t3762) * t3676 + mrSges(1,2);
t3632 = t3700 * t3682;
t3733 = pkin(1) * t3657 - pkin(6) * t3632;
t3597 = t3661 + t3632;
t3594 = pkin(1) + t3597;
t3765 = t3594 * t3830;
t3795 = t3677 * t3682;
t3768 = pkin(1) * t3795;
t3800 = t3676 * t3691;
t3736 = -(pkin(6) * t3800 + t3768) * t3833 + t3676 * t3765;
t3823 = ((-t3683 * t3718 + t3721 * t3692) * t3591 + (t3683 * t3721 + t3718 * t3692) * t3588) / (((t3621 * t3661 + t3733) * t3676 + t3575) * t3690 + t3736);
t3636 = t3700 * t3694;
t3601 = -t3685 * pkin(2) + t3636;
t3576 = t3601 * t3837;
t3623 = t3658 - pkin(6);
t3586 = -mrSges(3,1) * t3693 + mrSges(3,2) * t3684 + t3824;
t3717 = t3586 * t3694 - t3638 * t3685 - t3616;
t3720 = (-t3685 * t3586 - t3638 * t3694) * t3677 - (mrSges(3,1) * t3684 + mrSges(3,2) * t3693 + t3762) * t3676 + mrSges(1,2);
t3633 = t3700 * t3685;
t3732 = pkin(1) * t3658 - pkin(6) * t3633;
t3598 = t3663 + t3633;
t3595 = pkin(1) + t3598;
t3764 = t3595 * t3829;
t3793 = t3677 * t3685;
t3767 = pkin(1) * t3793;
t3799 = t3676 * t3694;
t3735 = -(pkin(6) * t3799 + t3767) * t3832 + t3676 * t3764;
t3822 = ((-t3686 * t3717 + t3720 * t3695) * t3592 + (t3686 * t3720 + t3717 * t3695) * t3589) / (((t3623 * t3663 + t3732) * t3676 + t3576) * t3693 + t3735);
t3637 = t3700 * t3697;
t3602 = -t3688 * pkin(2) + t3637;
t3577 = t3602 * t3837;
t3625 = t3659 - pkin(6);
t3587 = -mrSges(3,1) * t3696 + mrSges(3,2) * t3687 + t3824;
t3716 = t3587 * t3697 - t3638 * t3688 - t3616;
t3719 = (-t3688 * t3587 - t3638 * t3697) * t3677 - (mrSges(3,1) * t3687 + mrSges(3,2) * t3696 + t3762) * t3676 + mrSges(1,2);
t3634 = t3700 * t3688;
t3731 = pkin(1) * t3659 - pkin(6) * t3634;
t3599 = t3665 + t3634;
t3596 = pkin(1) + t3599;
t3763 = t3596 * t3828;
t3791 = t3677 * t3688;
t3766 = pkin(1) * t3791;
t3798 = t3676 * t3697;
t3734 = -(pkin(6) * t3798 + t3766) * t3831 + t3676 * t3763;
t3821 = ((-t3689 * t3716 + t3719 * t3698) * t3593 + (t3689 * t3719 + t3716 * t3698) * t3590) / (((t3625 * t3665 + t3731) * t3676 + t3577) * t3696 + t3734);
t3509 = -t3709 * t3638 + (-t3539 * t3682 + t3730 * t3691) * t3585;
t3627 = t3661 + pkin(1);
t3724 = -t3597 * pkin(6) + t3627 * t3657;
t3772 = pkin(6) * t3833;
t3790 = t3677 * t3690;
t3820 = t3509 / (pkin(1) * (-t3626 * t3682 + t3635) * t3790 + t3676 * (t3724 * t3690 - t3691 * t3772 + t3765));
t3510 = -t3708 * t3638 + (-t3540 * t3685 + t3729 * t3694) * t3586;
t3629 = t3663 + pkin(1);
t3723 = -t3598 * pkin(6) + t3629 * t3658;
t3771 = pkin(6) * t3832;
t3789 = t3677 * t3693;
t3819 = t3510 / (pkin(1) * (-t3628 * t3685 + t3636) * t3789 + t3676 * (t3723 * t3693 - t3694 * t3771 + t3764));
t3511 = -t3707 * t3638 + (-t3541 * t3688 + t3728 * t3697) * t3587;
t3631 = t3665 + pkin(1);
t3722 = -t3599 * pkin(6) + t3631 * t3659;
t3770 = pkin(6) * t3831;
t3788 = t3677 * t3696;
t3818 = t3511 / (pkin(1) * (-t3630 * t3688 + t3637) * t3788 + t3676 * (t3722 * t3696 - t3697 * t3770 + t3763));
t3817 = (pkin(1) + 0.2e1 * t3632) * t3626;
t3816 = (pkin(1) + 0.2e1 * t3633) * t3628;
t3815 = (pkin(1) + 0.2e1 * t3634) * t3630;
t3814 = (t3660 + t3668) * (t3660 + t3669);
t3813 = (t3662 + t3668) * (t3662 + t3669);
t3812 = (t3664 + t3668) * (t3664 + t3669);
t3811 = t3626 * t3676;
t3810 = t3628 * t3676;
t3809 = t3630 * t3676;
t3808 = t3668 * t3669;
t3807 = t3676 * t3677;
t3806 = t3676 * t3681;
t3805 = t3676 * t3682;
t3804 = t3676 * t3684;
t3803 = t3676 * t3685;
t3802 = t3676 * t3687;
t3801 = t3676 * t3688;
t3705 = 0.1e1 / pkin(3);
t3797 = t3676 * t3705;
t3796 = t3677 * t3681;
t3794 = t3677 * t3684;
t3792 = t3677 * t3687;
t3787 = t3677 * t3700;
t3786 = t3690 * t3691;
t3784 = t3693 * t3694;
t3782 = t3696 * t3697;
t3780 = mrSges(3,2) * t3834;
t3566 = t3692 * t3646 + t3683 * t3649;
t3779 = t3700 * t3566;
t3567 = -t3683 * t3646 + t3692 * t3649;
t3778 = t3700 * t3567;
t3568 = t3695 * t3647 + t3686 * t3650;
t3777 = t3700 * t3568;
t3569 = -t3686 * t3647 + t3695 * t3650;
t3776 = t3700 * t3569;
t3570 = t3698 * t3648 + t3689 * t3651;
t3775 = t3700 * t3570;
t3571 = -t3689 * t3648 + t3698 * t3651;
t3774 = t3700 * t3571;
t3773 = pkin(3) * t3837;
t3769 = -0.2e1 * (t3677 + 0.1e1) * (t3677 - 0.1e1);
t3615 = pkin(1) * t3787;
t3761 = t3681 * pkin(6) + pkin(3);
t3760 = t3684 * pkin(6) + pkin(3);
t3759 = t3687 * pkin(6) + pkin(3);
t3758 = t3626 * t3787;
t3757 = t3628 * t3787;
t3756 = t3630 * t3787;
t3755 = t3676 * t3787;
t3754 = t3676 * t3795;
t3753 = t3676 * t3793;
t3752 = t3676 * t3791;
t3751 = t3566 * t3806;
t3750 = t3567 * t3806;
t3749 = t3683 * t3814;
t3748 = t3568 * t3804;
t3747 = t3569 * t3804;
t3746 = t3686 * t3813;
t3745 = t3570 * t3802;
t3744 = t3571 * t3802;
t3743 = t3689 * t3812;
t3506 = ((t3727 * t3676 - t3834) * mrSges(3,1) + t3709 * mrSges(3,2)) * t3690 + (t3709 * mrSges(3,1) - t3727 * t3826 + t3780) * t3681;
t3584 = -t3676 * pkin(6) * t3700 - pkin(1) * t3836;
t3706 = pkin(2) ^ 2;
t3742 = t3506 / (t3615 * t3786 + (t3584 * t3690 - t3670 * t3773) * t3682 + ((pkin(2) * t3621 * t3690 - t3772) * t3691 + (pkin(1) * t3626 + pkin(2) * t3632 + t3706 * t3691) * t3681) * t3676) * t3797;
t3507 = ((t3726 * t3676 - t3834) * mrSges(3,1) + t3708 * mrSges(3,2)) * t3693 + (t3708 * mrSges(3,1) - t3726 * t3826 + t3780) * t3684;
t3741 = t3507 / (t3615 * t3784 + (t3584 * t3693 - t3672 * t3773) * t3685 + ((pkin(2) * t3623 * t3693 - t3771) * t3694 + (pkin(1) * t3628 + pkin(2) * t3633 + t3706 * t3694) * t3684) * t3676) * t3797;
t3508 = ((t3725 * t3676 - t3834) * mrSges(3,1) + t3707 * mrSges(3,2)) * t3696 + (t3707 * mrSges(3,1) - t3725 * t3826 + t3780) * t3687;
t3740 = t3508 / (t3615 * t3782 + (t3584 * t3696 - t3674 * t3773) * t3688 + ((pkin(2) * t3625 * t3696 - t3770) * t3697 + (pkin(1) * t3630 + pkin(2) * t3634 + t3706 * t3697) * t3687) * t3676) * t3797;
t3739 = t3692 * t3758;
t3738 = t3695 * t3757;
t3737 = t3698 * t3756;
t3715 = -(t3627 * t3682 + t3700 + (t3682 * t3660 - t3635) * t3691) * t3806 + (pkin(3) * t3786 + t3594) * t3790;
t3714 = -(t3629 * t3685 + t3700 + (t3685 * t3662 - t3636) * t3694) * t3804 + (pkin(3) * t3784 + t3595) * t3789;
t3713 = -(t3631 * t3688 + t3700 + (t3688 * t3664 - t3637) * t3697) * t3802 + (pkin(3) * t3782 + t3596) * t3788;
t3671 = t3691 ^ 2;
t3581 = t3671 * pkin(2) + t3682 * t3635 - pkin(2);
t3642 = t3671 - 0.2e1;
t3712 = ((t3642 * t3657 - pkin(6)) * t3690 + t3681 * t3581) * t3807 - (t3600 * t3690 + (t3761 - 0.2e1 * t3833) * t3682) * t3667 - t3682 * (-t3761 + t3833);
t3673 = t3694 ^ 2;
t3582 = t3673 * pkin(2) + t3685 * t3636 - pkin(2);
t3643 = t3673 - 0.2e1;
t3711 = ((t3643 * t3658 - pkin(6)) * t3693 + t3684 * t3582) * t3807 - (t3601 * t3693 + (t3760 - 0.2e1 * t3832) * t3685) * t3667 - t3685 * (-t3760 + t3832);
t3675 = t3697 ^ 2;
t3583 = t3675 * pkin(2) + t3688 * t3637 - pkin(2);
t3644 = t3675 - 0.2e1;
t3710 = ((t3644 * t3659 - pkin(6)) * t3696 + t3687 * t3583) * t3807 - (t3602 * t3696 + (t3759 - 0.2e1 * t3831) * t3688) * t3667 - t3688 * (-t3759 + t3831);
t3704 = pkin(3) ^ 2;
t3666 = pkin(1) * t3700;
t3624 = t3659 + pkin(6);
t3622 = t3658 + pkin(6);
t3620 = t3657 + pkin(6);
t3619 = pkin(1) * t3688 + t3700;
t3618 = pkin(1) * t3685 + t3700;
t3617 = pkin(1) * t3682 + t3700;
t3559 = -t3624 * t3801 + t3677 * t3630;
t3558 = -t3622 * t3803 + t3677 * t3628;
t3557 = -t3620 * t3805 + t3677 * t3626;
t3538 = 0.2e1 * t3689 * t3756 + t3698 * t3812;
t3537 = 0.2e1 * t3686 * t3757 + t3695 * t3813;
t3536 = 0.2e1 * t3683 * t3758 + t3692 * t3814;
t3535 = t3559 * t3698 + t3689 * t3619;
t3534 = t3558 * t3695 + t3686 * t3618;
t3533 = t3557 * t3692 + t3683 * t3617;
t3532 = -t3689 * t3559 + t3619 * t3698;
t3531 = -t3686 * t3558 + t3618 * t3695;
t3530 = -t3683 * t3557 + t3617 * t3692;
t3529 = -t3624 * t3809 + t3791 * t3812;
t3528 = -t3622 * t3810 + t3793 * t3813;
t3527 = -t3620 * t3811 + t3795 * t3814;
t3526 = t3529 * t3698 + t3689 * t3815;
t3525 = t3528 * t3695 + t3686 * t3816;
t3524 = t3527 * t3692 + t3683 * t3817;
t3523 = -t3529 * t3689 + t3698 * t3815;
t3522 = -t3528 * t3686 + t3695 * t3816;
t3521 = -t3527 * t3683 + t3692 * t3817;
t1 = [(-(-t3570 * t3791 + t3571 * t3697) * t3831 + (-pkin(3) * t3745 + (-pkin(2) * t3571 - t3677 * t3775) * t3697 + (t3570 * t3836 - t3774) * t3688) * t3696 - pkin(2) * t3745) * t3821 + (-t3710 * t3570 + t3713 * t3571) * t3818 - ((t3538 * t3651 + 0.2e1 * (t3737 - t3743 / 0.2e1) * t3648) * t3675 + (t3523 * t3651 - t3648 * t3526) * t3697 + (t3532 * t3651 - t3648 * t3535) * t3700) * t3740 + (-(-t3568 * t3793 + t3569 * t3694) * t3832 + (-pkin(3) * t3748 + (-pkin(2) * t3569 - t3677 * t3777) * t3694 + (t3568 * t3836 - t3776) * t3685) * t3693 - pkin(2) * t3748) * t3822 + (-t3711 * t3568 + t3714 * t3569) * t3819 - ((t3537 * t3650 + 0.2e1 * (t3738 - t3746 / 0.2e1) * t3647) * t3673 + (t3522 * t3650 - t3647 * t3525) * t3694 + (t3531 * t3650 - t3647 * t3534) * t3700) * t3741 + (-(-t3566 * t3795 + t3567 * t3691) * t3833 + (-pkin(3) * t3751 + (-pkin(2) * t3567 - t3677 * t3779) * t3691 + (t3566 * t3836 - t3778) * t3682) * t3690 - pkin(2) * t3751) * t3823 + (-t3712 * t3566 + t3715 * t3567) * t3820 - ((t3536 * t3649 + 0.2e1 * (t3739 - t3749 / 0.2e1) * t3646) * t3671 + (t3521 * t3649 - t3646 * t3524) * t3691 + (t3530 * t3649 - t3646 * t3533) * t3700) * t3742 - g(1) * m(4); (-(t3570 * t3697 + t3571 * t3791) * t3831 + (pkin(3) * t3744 + (-pkin(2) * t3570 + t3677 * t3774) * t3697 - (t3571 * t3836 + t3775) * t3688) * t3696 + pkin(2) * t3744) * t3821 + (t3713 * t3570 + t3710 * t3571) * t3818 - (((-0.2e1 * t3737 + t3743) * t3651 + t3538 * t3648) * t3675 + (t3523 * t3648 + t3526 * t3651) * t3697 + (t3532 * t3648 + t3535 * t3651) * t3700) * t3740 + (-(t3568 * t3694 + t3569 * t3793) * t3832 + (pkin(3) * t3747 + (-pkin(2) * t3568 + t3677 * t3776) * t3694 - (t3569 * t3836 + t3777) * t3685) * t3693 + pkin(2) * t3747) * t3822 + (t3714 * t3568 + t3711 * t3569) * t3819 - (((-0.2e1 * t3738 + t3746) * t3650 + t3537 * t3647) * t3673 + (t3522 * t3647 + t3525 * t3650) * t3694 + (t3531 * t3647 + t3534 * t3650) * t3700) * t3741 + (-(t3566 * t3691 + t3567 * t3795) * t3833 + (pkin(3) * t3750 + (-pkin(2) * t3566 + t3677 * t3778) * t3691 - (t3567 * t3836 + t3779) * t3682) * t3690 + pkin(2) * t3750) * t3823 + (t3715 * t3566 + t3712 * t3567) * t3820 - (((-0.2e1 * t3739 + t3749) * t3649 + t3536 * t3646) * t3671 + (t3521 * t3646 + t3524 * t3649) * t3691 + (t3530 * t3646 + t3533 * t3649) * t3700) * t3742 - g(2) * m(4); (-t3801 * t3831 + (-pkin(3) * t3792 + t3602 * t3676) * t3696 - pkin(2) * t3792) * t3821 + ((-t3602 * t3807 + t3827) * t3696 + (-pkin(6) * t3752 - t3583 * t3667 + t3697 * t3596) * t3687 + (-(t3644 * t3667 - t3675 + 0.1e1) * t3687 * t3696 + (0.2e1 * t3674 - 0.1e1) * t3752) * pkin(3)) / ((t3676 * t3722 + t3577) * t3696 + t3734) * t3511 + (-t3803 * t3832 + (-pkin(3) * t3794 + t3601 * t3676) * t3693 - pkin(2) * t3794) * t3822 + ((-t3601 * t3807 + t3827) * t3693 + (-pkin(6) * t3753 - t3582 * t3667 + t3694 * t3595) * t3684 + (-(t3643 * t3667 - t3673 + 0.1e1) * t3684 * t3693 + (0.2e1 * t3672 - 0.1e1) * t3753) * pkin(3)) / ((t3676 * t3723 + t3576) * t3693 + t3735) * t3510 + (-t3805 * t3833 + (-pkin(3) * t3796 + t3600 * t3676) * t3690 - pkin(2) * t3796) * t3823 + ((-t3600 * t3807 + t3827) * t3690 + (-pkin(6) * t3754 - t3581 * t3667 + t3691 * t3594) * t3681 + (-(t3642 * t3667 - t3671 + 0.1e1) * t3681 * t3690 + (0.2e1 * t3670 - 0.1e1) * t3754) * pkin(3)) / ((t3676 * t3724 + t3575) * t3690 + t3736) * t3509 - g(3) * m(4) + ((t3666 * t3697 + (-t3624 * t3677 * t3798 - t3619 + (t3675 * t3769 + t3667) * t3700) * t3630 + ((t3667 * t3812 - t3674 * t3704 + t3696 * t3825 - t3808) * t3697 - t3624 * t3755) * t3688) / ((t3615 * t3696 - (pkin(6) * t3696 - t3828) * t3809) * t3697 - t3696 * t3630 * t3766 + (t3731 * t3696 + (t3634 + pkin(1)) * t3828) * t3676) * t3508 + (t3666 * t3694 + (-t3622 * t3677 * t3799 - t3618 + (t3673 * t3769 + t3667) * t3700) * t3628 + ((t3667 * t3813 - t3672 * t3704 + t3693 * t3825 - t3808) * t3694 - t3622 * t3755) * t3685) / ((t3615 * t3693 - (pkin(6) * t3693 - t3829) * t3810) * t3694 - t3693 * t3628 * t3767 + (t3732 * t3693 + (t3633 + pkin(1)) * t3829) * t3676) * t3507 + (t3666 * t3691 + (-t3620 * t3677 * t3800 - t3617 + (t3671 * t3769 + t3667) * t3700) * t3626 + ((t3667 * t3814 - t3670 * t3704 + t3690 * t3825 - t3808) * t3691 - t3620 * t3755) * t3682) / ((t3690 * t3615 - (pkin(6) * t3690 - t3830) * t3811) * t3691 - t3690 * t3626 * t3768 + (t3733 * t3690 + (t3632 + pkin(1)) * t3830) * t3676) * t3506) * t3705;];
taugX  = t1;
