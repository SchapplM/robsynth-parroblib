% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR8V2G1A0
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
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR8V2G1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-07 13:12
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V2G1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:12:55
% EndTime: 2022-11-07 13:12:57
% DurationCPUTime: 1.95s
% Computational Cost: add. (918->188), mult. (1185->311), div. (123->6), fcn. (1119->33), ass. (0->141)
t3863 = MDP(3) - MDP(13);
t3806 = qJ(3,1) + pkin(5);
t3805 = qJ(3,2) + pkin(5);
t3804 = qJ(3,3) + pkin(5);
t3798 = qJ(2,3) + pkin(7);
t3780 = cos(t3798);
t3859 = pkin(3) * t3780;
t3799 = qJ(2,2) + pkin(7);
t3781 = cos(t3799);
t3858 = pkin(3) * t3781;
t3800 = qJ(2,1) + pkin(7);
t3782 = cos(t3800);
t3857 = pkin(3) * t3782;
t3856 = 0.2e1 * pkin(2) * pkin(3);
t3855 = 2 * pkin(1);
t3801 = legFrame(3,3);
t3783 = sin(t3801);
t3786 = cos(t3801);
t3808 = sin(qJ(1,3));
t3814 = cos(qJ(1,3));
t3746 = t3783 * t3814 + t3786 * t3808;
t3824 = t3783 * t3808 - t3786 * t3814;
t3720 = t3746 * g(1) + t3824 * g(2);
t3795 = -pkin(6) - t3804;
t3789 = 0.1e1 / t3795;
t3854 = t3720 * t3789;
t3810 = sin(qJ(1,2));
t3816 = cos(qJ(1,2));
t3764 = t3810 * g(1) - t3816 * g(2);
t3765 = t3816 * g(1) + t3810 * g(2);
t3802 = legFrame(2,3);
t3784 = sin(t3802);
t3787 = cos(t3802);
t3736 = t3764 * t3787 + t3765 * t3784;
t3796 = -pkin(6) - t3805;
t3790 = 0.1e1 / t3796;
t3853 = t3736 * t3790;
t3852 = t3746 * t3789;
t3851 = t3824 * t3789;
t3748 = t3784 * t3816 + t3787 * t3810;
t3850 = t3748 * t3790;
t3749 = -t3784 * t3810 + t3787 * t3816;
t3849 = t3749 * t3790;
t3803 = legFrame(1,3);
t3785 = sin(t3803);
t3788 = cos(t3803);
t3812 = sin(qJ(1,1));
t3818 = cos(qJ(1,1));
t3750 = t3785 * t3818 + t3788 * t3812;
t3797 = -pkin(6) - t3806;
t3791 = 0.1e1 / t3797;
t3848 = t3750 * t3791;
t3751 = -t3785 * t3812 + t3788 * t3818;
t3847 = t3751 * t3791;
t3777 = sin(t3798);
t3807 = sin(qJ(2,3));
t3761 = t3807 * pkin(2) + pkin(3) * t3777;
t3846 = t3761 * t3789;
t3778 = sin(t3799);
t3809 = sin(qJ(2,2));
t3762 = t3809 * pkin(2) + pkin(3) * t3778;
t3845 = t3762 * t3790;
t3779 = sin(t3800);
t3811 = sin(qJ(2,1));
t3763 = t3811 * pkin(2) + pkin(3) * t3779;
t3844 = t3763 * t3791;
t3843 = t3777 * t3854;
t3842 = t3780 * t3854;
t3721 = -t3750 * g(1) + t3751 * g(2);
t3841 = t3721 * t3779 * t3791;
t3722 = -t3748 * g(1) + t3749 * g(2);
t3840 = t3722 * t3778 * t3790;
t3752 = t3783 * g(1) - t3786 * g(2);
t3755 = t3786 * g(1) + t3783 * g(2);
t3726 = t3752 * t3814 + t3755 * t3808;
t3839 = t3726 * t3846;
t3754 = t3785 * g(1) - t3788 * g(2);
t3757 = t3788 * g(1) + t3785 * g(2);
t3728 = t3754 * t3818 + t3757 * t3812;
t3838 = t3728 * t3844;
t3837 = t3781 * t3853;
t3815 = cos(qJ(2,2));
t3836 = t3815 * t3853;
t3766 = g(1) * t3812 - g(2) * t3818;
t3767 = g(1) * t3818 + g(2) * t3812;
t3738 = t3766 * t3788 + t3767 * t3785;
t3835 = t3738 * t3782 * t3791;
t3834 = t3726 * t3852;
t3833 = t3726 * t3851;
t3753 = t3784 * g(1) - t3787 * g(2);
t3756 = t3787 * g(1) + t3784 * g(2);
t3727 = t3753 * t3816 + t3756 * t3810;
t3832 = t3727 * t3850;
t3831 = t3727 * t3849;
t3830 = t3728 * t3848;
t3829 = t3728 * t3847;
t3813 = cos(qJ(2,3));
t3792 = t3813 * pkin(2);
t3758 = 0.1e1 / (t3792 + t3859);
t3828 = t3758 * t3846;
t3793 = t3815 * pkin(2);
t3759 = 0.1e1 / (t3793 + t3858);
t3827 = t3759 * t3845;
t3817 = cos(qJ(2,1));
t3794 = t3817 * pkin(2);
t3760 = 0.1e1 / (t3794 + t3857);
t3826 = t3760 * t3844;
t3825 = t3720 * t3846 + g(3);
t3723 = t3808 * t3752 - t3755 * t3814;
t3724 = t3810 * t3753 - t3756 * t3816;
t3725 = t3812 * t3754 - t3757 * t3818;
t3823 = pkin(2) ^ 2;
t3822 = pkin(3) ^ 2;
t3821 = 0.2e1 * qJ(2,1);
t3820 = 0.2e1 * qJ(2,2);
t3819 = 0.2e1 * qJ(2,3);
t3776 = t3794 + pkin(1);
t3775 = t3793 + pkin(1);
t3774 = t3792 + pkin(1);
t3773 = t3775 * t3816;
t3772 = t3774 * t3814;
t3771 = t3818 * t3776;
t3770 = t3812 * t3776;
t3769 = t3810 * t3775;
t3768 = t3808 * t3774;
t3745 = -t3810 * t3796 + t3773;
t3744 = -t3808 * t3795 + t3772;
t3743 = -t3812 * t3797 + t3771;
t3742 = t3818 * t3797 + t3770;
t3741 = t3816 * t3796 + t3769;
t3740 = t3814 * t3795 + t3768;
t3739 = -t3766 * t3785 + t3767 * t3788;
t3737 = -t3764 * t3784 + t3765 * t3787;
t3735 = (t3814 * g(1) + t3808 * g(2)) * t3786 - (t3808 * g(1) - t3814 * g(2)) * t3783;
t3719 = -g(3) * t3817 + t3739 * t3811;
t3718 = -g(3) * t3815 + t3737 * t3809;
t3717 = -g(3) * t3813 + t3735 * t3807;
t3716 = (-t3806 * t3818 + t3770) * t3757 + t3754 * (t3806 * t3812 + t3771);
t3715 = (-t3805 * t3816 + t3769) * t3756 + t3753 * (t3805 * t3810 + t3773);
t3714 = (-t3804 * t3814 + t3768) * t3755 + t3752 * (t3804 * t3808 + t3772);
t1 = [(-t3829 - t3831 + t3833) * MDP(2) + (-t3749 * t3836 + t3813 * t3833 - t3817 * t3829) * MDP(9) + (-t3807 * t3833 + t3809 * t3831 + t3811 * t3829) * MDP(10) + (-t3749 * t3837 - t3751 * t3835 + t3824 * t3842) * MDP(11) + (-t3749 * t3840 - t3751 * t3841 - t3824 * t3843) * MDP(12) + (-(t3751 * t3716 - (-t3785 * t3742 + t3743 * t3788 + t3751 * t3857) * t3728) * t3791 - (t3749 * t3715 - (-t3784 * t3741 + t3745 * t3787 + t3749 * t3858) * t3727) * t3790 - (-t3824 * t3714 - (-t3783 * t3740 + t3744 * t3786 - t3824 * t3859) * t3726) * t3789) * MDP(14) - g(1) * MDP(15) - t3863 * (t3723 * t3851 - t3724 * t3849 - t3725 * t3847); (-t3830 - t3832 - t3834) * MDP(2) + (-t3748 * t3836 - t3813 * t3834 - t3817 * t3830) * MDP(9) + (t3807 * t3834 + t3809 * t3832 + t3811 * t3830) * MDP(10) + (-t3746 * t3842 - t3748 * t3837 - t3750 * t3835) * MDP(11) + (t3746 * t3843 - t3748 * t3840 - t3750 * t3841) * MDP(12) + (-(t3750 * t3716 - (t3742 * t3788 + t3743 * t3785 + t3750 * t3857) * t3728) * t3791 - (t3748 * t3715 - (t3741 * t3787 + t3745 * t3784 + t3748 * t3858) * t3727) * t3790 - (t3746 * t3714 - (t3740 * t3786 + t3744 * t3783 + t3746 * t3859) * t3726) * t3789) * MDP(14) - g(2) * MDP(15) + t3863 * (t3723 * t3852 + t3724 * t3850 + t3725 * t3848); (-t3726 * t3828 - t3727 * t3827 - t3728 * t3826) * MDP(2) + ((-t3817 * t3838 + t3719) * t3760 + (-t3762 * t3836 + t3718) * t3759 + (-t3813 * t3839 + t3717) * t3758) * MDP(9) + ((t3739 * t3817 + (g(3) + t3838) * t3811) * t3760 + (t3737 * t3815 + (t3727 * t3845 + g(3)) * t3809) * t3759 + (t3735 * t3813 + (g(3) + t3839) * t3807) * t3758) * MDP(10) + ((t3739 * t3779 + (-t3738 * t3844 - g(3)) * t3782) * t3760 + (t3737 * t3778 + (-t3736 * t3845 - g(3)) * t3781) * t3759 + (t3735 * t3777 - t3825 * t3780) * t3758) * MDP(11) + ((t3739 * t3782 + (-t3721 * t3844 + g(3)) * t3779) * t3760 + (t3737 * t3781 + (-t3722 * t3845 + g(3)) * t3778) * t3759 + (t3735 * t3780 + t3825 * t3777) * t3758) * MDP(12) + (-t3716 * t3826 + (sin(t3821 + pkin(7)) * t3856 + t3822 * sin(0.2e1 * t3800) + t3823 * sin(t3821) + t3763 * t3855) * t3791 * t3760 * t3728 / 0.2e1 - t3715 * t3827 + (sin(t3820 + pkin(7)) * t3856 + t3822 * sin(0.2e1 * t3799) + t3823 * sin(t3820) + t3762 * t3855) * t3790 * t3759 * t3727 / 0.2e1 - t3714 * t3828 + (sin(t3819 + pkin(7)) * t3856 + t3822 * sin(0.2e1 * t3798) + t3823 * sin(t3819) + t3761 * t3855) * t3789 * t3758 * t3726 / 0.2e1 + (t3758 * t3717 + t3759 * t3718 + t3760 * t3719) * pkin(2)) * MDP(14) - g(3) * MDP(15) + t3863 * (t3723 * t3828 + t3724 * t3827 + t3725 * t3826);];
taugX  = t1;
