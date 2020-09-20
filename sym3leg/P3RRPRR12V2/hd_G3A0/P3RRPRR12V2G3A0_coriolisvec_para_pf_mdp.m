% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RRPRR12V2G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR12V2G3A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR12V2G3A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:30:00
% EndTime: 2020-08-06 19:30:35
% DurationCPUTime: 35.04s
% Computational Cost: add. (242964->863), mult. (327555->1640), div. (12789->12), fcn. (169938->18), ass. (0->594)
t3886 = xDP(2);
t3889 = pkin(2) + pkin(3);
t4150 = (t3886 * t3889);
t3818 = pkin(1) * t4150;
t3887 = xDP(1);
t4149 = (t3887 * t3889);
t3819 = pkin(1) * t4149;
t3870 = legFrame(3,2);
t3841 = sin(t3870);
t3844 = cos(t3870);
t3879 = cos(qJ(2,3));
t3862 = t3879 ^ 2;
t3873 = sin(qJ(2,3));
t3880 = cos(qJ(1,3));
t3885 = xDP(3);
t3888 = pkin(5) - pkin(6);
t3832 = t3888 * t3880;
t3874 = sin(qJ(1,3));
t3829 = t3874 * t3888;
t4117 = pkin(1) * t3880 + t3829;
t4326 = -t3841 * t3886 + t3844 * t3887;
t3914 = -(pkin(1) * t3874 - t3832) * t3885 + t4326 * t4117;
t3947 = t3841 * t3887 + t3844 * t3886;
t4178 = (qJ(3,3) + t3889) * (-qJ(3,3) + t3889);
t4030 = t3880 * t4178;
t4031 = t3874 * t4178;
t4102 = -2 * t4149;
t4103 = -2 * t4150;
t4154 = t3885 * t3874;
t4279 = qJ(3,3) * t3887;
t4280 = qJ(3,3) * t3886;
t3867 = pkin(1) * qJ(3,3);
t4032 = t3873 * t4178;
t4329 = t4032 - t3867;
t3838 = t3873 * qJ(3,3);
t4353 = 0.2e1 * t3838;
t3719 = ((qJ(3,3) * t4103 + t3887 * t4030) * t3844 + (qJ(3,3) * t4102 - t3886 * t4030) * t3841 - t3885 * t4031) * t3862 + (((t3880 * t4326 - t4154) * t4353 + t3914) * t3889 + t4329 * t3947) * t3879 + (qJ(3,3) * t3914 + t3818 * t3844 + t3819 * t3841) * t3873 - ((-t3880 * t4279 - t4150) * t3844 + (t3880 * t4280 - t4149) * t3841 + qJ(3,3) * t4154) * qJ(3,3);
t3821 = t3838 + pkin(1);
t4165 = t3879 * t3889;
t3803 = t4165 + t3821;
t3794 = 0.1e1 / t3803;
t3891 = 0.1e1 / qJ(3,3);
t4213 = t3794 * t3891;
t3713 = t3719 * t4213;
t4311 = pkin(1) * t3873;
t3826 = qJ(3,3) + t4311;
t3950 = t3821 * t3880 + t3829;
t4151 = t3885 * t3889;
t4334 = t3821 * t3874 - t3832;
t3740 = ((t3880 * t4149 - t4280) * t3844 + (-t3880 * t4150 - t4279) * t3841 - t3874 * t4151) * t3862 + ((t3873 * t4150 + t3887 * t3950) * t3844 + (t3873 * t4149 - t3886 * t3950) * t3841 - t4334 * t3885) * t3879 + t3947 * t3826;
t4064 = t3740 * t4213;
t3734 = pkin(2) * t4064;
t3704 = t3734 - t3713;
t4017 = pkin(3) * t4064;
t3689 = -t4017 - t3704;
t3773 = -t3874 * t4326 - t3880 * t3885;
t4227 = t3773 * t3888;
t4051 = t3873 * t4227;
t4147 = t3889 * t3891;
t3731 = (-t3740 * t4147 + t4051) * t3794;
t3795 = 0.1e1 / t3803 ^ 2;
t3890 = qJ(3,3) ^ 2;
t4114 = pkin(1) ^ 2 + pkin(5) ^ 2;
t3969 = t4114 + (-0.2e1 * pkin(5) + pkin(6)) * pkin(6);
t3809 = t3890 + t3969;
t3892 = 0.1e1 / qJ(3,3) ^ 2;
t4340 = t3821 * t3879;
t3941 = 0.2e1 * t4340;
t4229 = t3773 * t3794;
t4053 = t3862 * t4229;
t4009 = qJ(3,3) * t4053;
t4142 = 0.2e1 * t3867;
t4123 = t3888 * t3689 + t4142 * t4229;
t4164 = t3879 * t3891;
t4248 = t3740 * t3891;
t4249 = t3740 * t3888;
t3659 = ((t4123 * t3873 + (t3879 * t4249 + (t3862 * t4178 + t3889 * t3941 + t3809) * t3773) * t3794) * t3773 * t4164 + (-(-t3888 * t4009 + ((-pkin(3) * t4248 + t4051) * t3794 - t3704) * t4165 + t3689 * t3821) * t3740 - (-t3731 * t3879 + t3821 * t4064) * t3719) * t3892) * t3795;
t4308 = pkin(2) * t3659;
t4256 = t3689 * t3889;
t4284 = qJ(3,3) * t3740;
t3680 = -t3794 * t4284 + t4256;
t3835 = 0.2e1 * t3862 - 0.1e1;
t3866 = t3889 ^ 2;
t4214 = t3794 * t3873;
t4052 = t3773 * t4214;
t4113 = -0.4e1 * pkin(1) * t3889;
t4148 = t3888 * t3889;
t4347 = 0.2e1 * pkin(1);
t4348 = ((((t3794 * t4032 * t4227 + t3680 * t3889) * t3879 + pkin(1) * t4256 + (t3680 * t3873 + (-t3835 * t3773 * t4148 - pkin(1) * t3740) * t3794) * qJ(3,3)) * t3740 + (-t3731 * t4165 + ((pkin(1) * t3891 + t3873) * t3889 * t3740 + (t3862 - 0.1e1) * qJ(3,3) * t4227) * t3794) * t3719) * t3892 + (-(t3866 - 0.3e1 * t3890) * t4053 * t4165 + (t3888 * (t3713 - 0.2e1 * t3734 - 0.2e1 * t4017) * qJ(3,3) + (-0.3e1 * (-t3890 / 0.3e1 + t3866) * t3838 + (t3890 - t3866) * t4347) * t4229) * t3862 + (-t3689 * t3873 * t4148 + (-t3773 * (0.3e1 * t3890 + t3969) * t3889 + (t3773 * t4113 - t4249) * t3838) * t3794) * t3879 - qJ(3,3) * (t3809 * t4052 + t4123)) * t3773 * t3891) * t3795;
t4356 = t4348 + t4308;
t3871 = legFrame(2,2);
t3842 = sin(t3871);
t3845 = cos(t3871);
t3881 = cos(qJ(2,2));
t3863 = t3881 ^ 2;
t3875 = sin(qJ(2,2));
t3882 = cos(qJ(1,2));
t3833 = t3888 * t3882;
t3876 = sin(qJ(1,2));
t3830 = t3876 * t3888;
t4116 = pkin(1) * t3882 + t3830;
t4327 = -t3842 * t3886 + t3845 * t3887;
t3913 = -(pkin(1) * t3876 - t3833) * t3885 + t4327 * t4116;
t3945 = t3842 * t3887 + t3845 * t3886;
t4177 = (qJ(3,2) + t3889) * (-qJ(3,2) + t3889);
t4027 = t3882 * t4177;
t4028 = t3876 * t4177;
t4153 = t3885 * t3876;
t4285 = qJ(3,2) * t3887;
t4286 = qJ(3,2) * t3886;
t3868 = pkin(1) * qJ(3,2);
t4029 = t3875 * t4177;
t4330 = t4029 - t3868;
t3839 = t3875 * qJ(3,2);
t4352 = 0.2e1 * t3839;
t3720 = ((qJ(3,2) * t4103 + t3887 * t4027) * t3845 + (qJ(3,2) * t4102 - t3886 * t4027) * t3842 - t3885 * t4028) * t3863 + (((t3882 * t4327 - t4153) * t4352 + t3913) * t3889 + t4330 * t3945) * t3881 + (qJ(3,2) * t3913 + t3818 * t3845 + t3819 * t3842) * t3875 - ((-t3882 * t4285 - t4150) * t3845 + (t3882 * t4286 - t4149) * t3842 + qJ(3,2) * t4153) * qJ(3,2);
t3823 = t3839 + pkin(1);
t4161 = t3881 * t3889;
t3804 = t4161 + t3823;
t3797 = 0.1e1 / t3804;
t3894 = 0.1e1 / qJ(3,2);
t4208 = t3797 * t3894;
t3714 = t3720 * t4208;
t4310 = pkin(1) * t3875;
t3827 = qJ(3,2) + t4310;
t3949 = t3823 * t3882 + t3830;
t4333 = t3823 * t3876 - t3833;
t3741 = ((t3882 * t4149 - t4286) * t3845 + (-t3882 * t4150 - t4285) * t3842 - t3876 * t4151) * t3863 + ((t3875 * t4150 + t3887 * t3949) * t3845 + (t3875 * t4149 - t3886 * t3949) * t3842 - t4333 * t3885) * t3881 + t3945 * t3827;
t4061 = t3741 * t4208;
t3735 = pkin(2) * t4061;
t3706 = t3735 - t3714;
t4016 = pkin(3) * t4061;
t3690 = -t4016 - t3706;
t3774 = -t3876 * t4327 - t3882 * t3885;
t4224 = t3774 * t3888;
t4048 = t3875 * t4224;
t4146 = t3889 * t3894;
t3732 = (-t3741 * t4146 + t4048) * t3797;
t3798 = 0.1e1 / t3804 ^ 2;
t3893 = qJ(3,2) ^ 2;
t3810 = t3893 + t3969;
t3895 = 0.1e1 / qJ(3,2) ^ 2;
t4339 = t3823 * t3881;
t3940 = 0.2e1 * t4339;
t4226 = t3774 * t3797;
t4050 = t3863 * t4226;
t4010 = qJ(3,2) * t4050;
t4143 = 0.2e1 * t3868;
t4122 = t3888 * t3690 + t4143 * t4226;
t4160 = t3881 * t3894;
t4245 = t3741 * t3894;
t4246 = t3741 * t3888;
t3660 = ((t4122 * t3875 + (t3881 * t4246 + (t3863 * t4177 + t3889 * t3940 + t3810) * t3774) * t3797) * t3774 * t4160 + (-(-t3888 * t4010 + ((-pkin(3) * t4245 + t4048) * t3797 - t3706) * t4161 + t3690 * t3823) * t3741 - (-t3732 * t3881 + t3823 * t4061) * t3720) * t3895) * t3798;
t4307 = pkin(2) * t3660;
t4255 = t3690 * t3889;
t4290 = qJ(3,2) * t3741;
t3681 = -t3797 * t4290 + t4255;
t3836 = 0.2e1 * t3863 - 0.1e1;
t4209 = t3797 * t3875;
t4049 = t3774 * t4209;
t4349 = ((((t3797 * t4029 * t4224 + t3681 * t3889) * t3881 + pkin(1) * t4255 + (t3681 * t3875 + (-t3836 * t3774 * t4148 - pkin(1) * t3741) * t3797) * qJ(3,2)) * t3741 + (-t3732 * t4161 + ((pkin(1) * t3894 + t3875) * t3889 * t3741 + (t3863 - 0.1e1) * qJ(3,2) * t4224) * t3797) * t3720) * t3895 + (-(t3866 - 0.3e1 * t3893) * t4050 * t4161 + (t3888 * (t3714 - 0.2e1 * t3735 - 0.2e1 * t4016) * qJ(3,2) + (-0.3e1 * (-t3893 / 0.3e1 + t3866) * t3839 + (t3893 - t3866) * t4347) * t4226) * t3863 + (-t3690 * t3875 * t4148 + (-t3774 * (0.3e1 * t3893 + t3969) * t3889 + (t3774 * t4113 - t4246) * t3839) * t3797) * t3881 - qJ(3,2) * (t3810 * t4049 + t4122)) * t3774 * t3894) * t3798;
t4355 = t4349 + t4307;
t3872 = legFrame(1,2);
t3843 = sin(t3872);
t3846 = cos(t3872);
t3883 = cos(qJ(2,1));
t3864 = t3883 ^ 2;
t3877 = sin(qJ(2,1));
t3884 = cos(qJ(1,1));
t3834 = t3888 * t3884;
t3878 = sin(qJ(1,1));
t3831 = t3878 * t3888;
t4115 = pkin(1) * t3884 + t3831;
t4328 = -t3843 * t3886 + t3846 * t3887;
t3912 = -(pkin(1) * t3878 - t3834) * t3885 + t4328 * t4115;
t3943 = t3843 * t3887 + t3846 * t3886;
t4176 = (qJ(3,1) + t3889) * (-qJ(3,1) + t3889);
t4024 = t3884 * t4176;
t4025 = t3878 * t4176;
t4152 = t3885 * t3878;
t4291 = qJ(3,1) * t3887;
t4292 = qJ(3,1) * t3886;
t3869 = pkin(1) * qJ(3,1);
t4026 = t3877 * t4176;
t4331 = t4026 - t3869;
t3840 = t3877 * qJ(3,1);
t4351 = 0.2e1 * t3840;
t3721 = ((qJ(3,1) * t4103 + t3887 * t4024) * t3846 + (qJ(3,1) * t4102 - t3886 * t4024) * t3843 - t3885 * t4025) * t3864 + (((t3884 * t4328 - t4152) * t4351 + t3912) * t3889 + t4331 * t3943) * t3883 + (qJ(3,1) * t3912 + t3818 * t3846 + t3819 * t3843) * t3877 - ((-t3884 * t4291 - t4150) * t3846 + (t3884 * t4292 - t4149) * t3843 + qJ(3,1) * t4152) * qJ(3,1);
t3825 = t3840 + pkin(1);
t4157 = t3883 * t3889;
t3805 = t4157 + t3825;
t3800 = 0.1e1 / t3805;
t3897 = 0.1e1 / qJ(3,1);
t4203 = t3800 * t3897;
t3715 = t3721 * t4203;
t4309 = pkin(1) * t3877;
t3828 = qJ(3,1) + t4309;
t3948 = t3825 * t3884 + t3831;
t4332 = t3825 * t3878 - t3834;
t3742 = ((t3884 * t4149 - t4292) * t3846 + (-t3884 * t4150 - t4291) * t3843 - t3878 * t4151) * t3864 + ((t3877 * t4150 + t3887 * t3948) * t3846 + (t3877 * t4149 - t3886 * t3948) * t3843 - t4332 * t3885) * t3883 + t3943 * t3828;
t4058 = t3742 * t4203;
t3736 = pkin(2) * t4058;
t3708 = t3736 - t3715;
t4015 = pkin(3) * t4058;
t3691 = -t4015 - t3708;
t3775 = -t3878 * t4328 - t3884 * t3885;
t4221 = t3775 * t3888;
t4045 = t3877 * t4221;
t4145 = t3889 * t3897;
t3733 = (-t3742 * t4145 + t4045) * t3800;
t3801 = 0.1e1 / t3805 ^ 2;
t3896 = qJ(3,1) ^ 2;
t3811 = t3896 + t3969;
t3898 = 0.1e1 / qJ(3,1) ^ 2;
t4338 = t3825 * t3883;
t3939 = 0.2e1 * t4338;
t4223 = t3775 * t3800;
t4047 = t3864 * t4223;
t4011 = qJ(3,1) * t4047;
t4144 = 0.2e1 * t3869;
t4121 = t3888 * t3691 + t4144 * t4223;
t4156 = t3883 * t3897;
t4242 = t3742 * t3897;
t4243 = t3742 * t3888;
t3661 = ((t4121 * t3877 + (t3883 * t4243 + (t3864 * t4176 + t3889 * t3939 + t3811) * t3775) * t3800) * t3775 * t4156 + (-(-t3888 * t4011 + ((-pkin(3) * t4242 + t4045) * t3800 - t3708) * t4157 + t3691 * t3825) * t3742 - (-t3733 * t3883 + t3825 * t4058) * t3721) * t3898) * t3801;
t4306 = pkin(2) * t3661;
t4254 = t3691 * t3889;
t4296 = qJ(3,1) * t3742;
t3682 = -t3800 * t4296 + t4254;
t3837 = 0.2e1 * t3864 - 0.1e1;
t4204 = t3800 * t3877;
t4046 = t3775 * t4204;
t4350 = ((((t3800 * t4026 * t4221 + t3682 * t3889) * t3883 + pkin(1) * t4254 + (t3682 * t3877 + (-t3837 * t3775 * t4148 - pkin(1) * t3742) * t3800) * qJ(3,1)) * t3742 + (-t3733 * t4157 + ((pkin(1) * t3897 + t3877) * t3889 * t3742 + (t3864 - 0.1e1) * qJ(3,1) * t4221) * t3800) * t3721) * t3898 + (-(t3866 - 0.3e1 * t3896) * t4047 * t4157 + (t3888 * (t3715 - 0.2e1 * t3736 - 0.2e1 * t4015) * qJ(3,1) + (-0.3e1 * (-t3896 / 0.3e1 + t3866) * t3840 + (t3896 - t3866) * t4347) * t4223) * t3864 + (-t3691 * t3877 * t4148 + (-t3775 * (0.3e1 * t3896 + t3969) * t3889 + (t3775 * t4113 - t4243) * t3840) * t3800) * t3883 - qJ(3,1) * (t3811 * t4046 + t4121)) * t3775 * t3897) * t3801;
t4354 = t4350 + t4306;
t4346 = 0.2e1 * t3841;
t4345 = 0.2e1 * t3842;
t4344 = 0.2e1 * t3843;
t4212 = t3795 * t3891;
t4063 = t3740 * t4212;
t4020 = 0.2e1 * t4063;
t4343 = t3794 * t4020;
t4207 = t3798 * t3894;
t4060 = t3741 * t4207;
t4019 = 0.2e1 * t4060;
t4342 = t3797 * t4019;
t4202 = t3801 * t3897;
t4057 = t3742 * t4202;
t4018 = 0.2e1 * t4057;
t4341 = t3800 * t4018;
t4300 = pkin(2) * t3883;
t4337 = (pkin(1) + t4300) * t3877;
t4301 = pkin(2) * t3881;
t4336 = (pkin(1) + t4301) * t3875;
t4302 = pkin(2) * t3879;
t4335 = (pkin(1) + t4302) * t3873;
t4321 = -0.2e1 * t3704;
t4320 = -0.2e1 * t3706;
t4319 = -0.2e1 * t3708;
t4318 = -0.2e1 * t3844;
t4317 = -0.2e1 * t3845;
t4316 = -0.2e1 * t3846;
t4315 = 0.4e1 * t3862;
t4314 = 0.4e1 * t3863;
t4313 = 0.4e1 * t3864;
t4312 = -0.2e1 * t3889;
t4305 = pkin(2) * t3862;
t4304 = pkin(2) * t3863;
t4303 = pkin(2) * t3864;
t3796 = t3794 * t3795;
t4210 = t3796 * t3891;
t4250 = t3740 * t3773;
t3984 = t4210 * t4250;
t4173 = t3873 * t3891;
t4040 = t3796 * t4173;
t4174 = t3873 * t3889;
t4228 = t3773 * t3795;
t4281 = qJ(3,3) * t3879;
t3674 = -(t3713 * t3873 + (t4227 + (-t3873 * t4147 + t3879) * t3740) * t3794) * t4228 + (t4174 - t4281) * t3984 - t3773 * t3719 * t4040;
t4299 = pkin(5) * t3674;
t3799 = t3797 * t3798;
t4205 = t3799 * t3894;
t4247 = t3741 * t3774;
t3983 = t4205 * t4247;
t4170 = t3875 * t3894;
t4037 = t3799 * t4170;
t4171 = t3875 * t3889;
t4225 = t3774 * t3798;
t4287 = qJ(3,2) * t3881;
t3675 = -(t3714 * t3875 + (t4224 + (-t3875 * t4146 + t3881) * t3741) * t3797) * t4225 + (t4171 - t4287) * t3983 - t3774 * t3720 * t4037;
t4298 = pkin(5) * t3675;
t3802 = t3800 * t3801;
t4200 = t3802 * t3897;
t4244 = t3742 * t3775;
t3982 = t4200 * t4244;
t4167 = t3877 * t3897;
t4034 = t3802 * t4167;
t4168 = t3877 * t3889;
t4222 = t3775 * t3801;
t4293 = qJ(3,1) * t3883;
t3676 = -(t3715 * t3877 + (t4221 + (-t3877 * t4145 + t3883) * t3742) * t3800) * t4222 + (t4168 - t4293) * t3982 - t3775 * t3721 * t4034;
t4297 = pkin(5) * t3676;
t4295 = qJ(3,1) * t3843;
t4294 = qJ(3,1) * t3846;
t4289 = qJ(3,2) * t3842;
t4288 = qJ(3,2) * t3845;
t4283 = qJ(3,3) * t3841;
t4282 = qJ(3,3) * t3844;
t4278 = t3862 * qJ(3,3);
t4277 = t3863 * qJ(3,2);
t4276 = t3864 * qJ(3,1);
t4274 = t3659 * t3879;
t4273 = t3660 * t3881;
t4272 = t3661 * t3883;
t4271 = t3674 * t3794;
t4270 = t3674 * t3873;
t4269 = t3674 * t3874;
t4268 = t3674 * t3879;
t4267 = t3674 * t3891;
t4266 = t3675 * t3797;
t4265 = t3675 * t3875;
t4264 = t3675 * t3876;
t4263 = t3675 * t3881;
t4262 = t3675 * t3894;
t4261 = t3676 * t3800;
t4260 = t3676 * t3877;
t4259 = t3676 * t3878;
t4258 = t3676 * t3883;
t4257 = t3676 * t3897;
t3737 = t3740 ^ 2;
t4253 = t3737 * t3892;
t3738 = t3741 ^ 2;
t4252 = t3738 * t3895;
t3739 = t3742 ^ 2;
t4251 = t3739 * t3898;
t4084 = t3880 * t3838;
t3783 = t4084 + t4117;
t4163 = t3880 * t3889;
t3752 = (-t3841 * t4163 - t4282) * t3862 + (-t3783 * t3841 + t3844 * t4174) * t3879 + t3844 * t3826;
t4241 = t3752 * t3891;
t3753 = (t3844 * t4163 - t4283) * t3862 + (t3783 * t3844 + t3841 * t4174) * t3879 + t3841 * t3826;
t4240 = t3753 * t3891;
t4085 = t3882 * t3839;
t3785 = t4085 + t4116;
t4159 = t3882 * t3889;
t3754 = (-t3842 * t4159 - t4288) * t3863 + (-t3785 * t3842 + t3845 * t4171) * t3881 + t3845 * t3827;
t4239 = t3754 * t3894;
t3755 = (t3845 * t4159 - t4289) * t3863 + (t3785 * t3845 + t3842 * t4171) * t3881 + t3842 * t3827;
t4238 = t3755 * t3894;
t4086 = t3884 * t3840;
t3787 = t4086 + t4115;
t4155 = t3884 * t3889;
t3756 = (-t3843 * t4155 - t4294) * t3864 + (-t3787 * t3843 + t3846 * t4168) * t3883 + t3846 * t3828;
t4237 = t3756 * t3897;
t3757 = (t3846 * t4155 - t4295) * t3864 + (t3787 * t3846 + t3843 * t4168) * t3883 + t3843 * t3828;
t4236 = t3757 * t3897;
t3820 = t4353 + pkin(1);
t3758 = -t3862 * t4031 - (t3820 * t3874 - t3832) * t4165 - qJ(3,3) * (t3826 * t3874 - t3873 * t3832);
t4235 = t3758 * t3873;
t3822 = t4352 + pkin(1);
t3759 = -t3863 * t4028 - (t3822 * t3876 - t3833) * t4161 - qJ(3,2) * (t3827 * t3876 - t3875 * t3833);
t4234 = t3759 * t3875;
t3824 = t4351 + pkin(1);
t3760 = -t3864 * t4025 - (t3824 * t3878 - t3834) * t4157 - qJ(3,1) * (t3828 * t3878 - t3877 * t3834);
t4233 = t3760 * t3877;
t3770 = t3773 ^ 2;
t3767 = t3770 * t3795;
t4232 = t3770 * t3862;
t3771 = t3774 ^ 2;
t3768 = t3771 * t3798;
t4231 = t3771 * t3863;
t3772 = t3775 ^ 2;
t3769 = t3772 * t3801;
t4230 = t3772 * t3864;
t3776 = t3874 * t4165 + t4334;
t4220 = t3776 * t3879;
t4219 = t3776 * t3891;
t3777 = t3876 * t4161 + t4333;
t4218 = t3777 * t3881;
t4217 = t3777 * t3894;
t3778 = t3878 * t4157 + t4332;
t4216 = t3778 * t3883;
t4215 = t3778 * t3897;
t4211 = t3795 * t3892;
t4206 = t3798 * t3895;
t4201 = t3801 * t3898;
t4199 = t3820 * t3879;
t4198 = t3822 * t3881;
t4197 = t3824 * t3883;
t4196 = t3841 * t3874;
t4194 = t3841 * t3889;
t4193 = t3842 * t3876;
t4191 = t3842 * t3889;
t4190 = t3843 * t3878;
t4188 = t3843 * t3889;
t4187 = t3844 * t3874;
t4185 = t3844 * t3889;
t4184 = t3845 * t3876;
t4182 = t3845 * t3889;
t4181 = t3846 * t3878;
t4179 = t3846 * t3889;
t4175 = t3873 * t3879;
t4172 = t3875 * t3881;
t4169 = t3877 * t3883;
t4166 = t3879 * t3880;
t4162 = t3881 * t3882;
t4158 = t3883 * t3884;
t3701 = t3734 - t3713 / 0.2e1;
t3705 = t3734 - 0.2e1 * t3713;
t3901 = pkin(2) ^ 2;
t3847 = -t3890 + t3901;
t3936 = qJ(3,3) * t4274 - t4356 * t3873;
t4023 = t4229 * t4321;
t4041 = t3794 * t4164;
t4141 = t3936 * pkin(5) + (pkin(2) * t3941 + t3847 * t3862 + t3873 * t4142 + t3890 + t4114) * t3674 + 0.4e1 * t3701 * t4009 + (-t3740 * pkin(5) * t3705 + 0.2e1 * (-(-pkin(2) * t3719 + t3740 * t3847) * t3873 + pkin(1) * t4284) * t4229) * t4041 + (-pkin(5) * t3737 * t4212 + pkin(1) * t4023) * t3873 + qJ(3,3) * t4023;
t3702 = t3735 - t3714 / 0.2e1;
t3707 = t3735 - 0.2e1 * t3714;
t3848 = -t3893 + t3901;
t3937 = qJ(3,2) * t4273 - t4355 * t3875;
t4022 = t4226 * t4320;
t4038 = t3797 * t4160;
t4140 = t3937 * pkin(5) + (pkin(2) * t3940 + t3848 * t3863 + t3875 * t4143 + t3893 + t4114) * t3675 + 0.4e1 * t3702 * t4010 + (-t3741 * pkin(5) * t3707 + 0.2e1 * (-(-pkin(2) * t3720 + t3741 * t3848) * t3875 + pkin(1) * t4290) * t4226) * t4038 + (-pkin(5) * t3738 * t4207 + pkin(1) * t4022) * t3875 + qJ(3,2) * t4022;
t3703 = t3736 - t3715 / 0.2e1;
t3709 = t3736 - 0.2e1 * t3715;
t3849 = -t3896 + t3901;
t3938 = qJ(3,1) * t4272 - t4354 * t3877;
t4021 = t4223 * t4319;
t4035 = t3800 * t4156;
t4139 = t3938 * pkin(5) + (pkin(2) * t3939 + t3849 * t3864 + t3877 * t4144 + t3896 + t4114) * t3676 + 0.4e1 * t3703 * t4011 + (-t3742 * pkin(5) * t3709 + 0.2e1 * (-(-pkin(2) * t3721 + t3742 * t3849) * t3877 + pkin(1) * t4296) * t4223) * t4035 + (-pkin(5) * t3739 * t4202 + pkin(1) * t4021) * t3877 + qJ(3,1) * t4021;
t3812 = -pkin(2) * t3873 + t4281;
t4109 = -0.2e1 * t4278;
t4120 = 0.2e1 * t3740 * t3719 * t4211 + pkin(2) * t3767;
t4138 = (t3901 + t3890) * t3659 + t3812 * t4299 + pkin(2) * t4348 + t4120 * qJ(3,3) + (-(-t3847 * t3873 + t3867) * t3879 + (t4109 + t4311) * pkin(2)) * t3767;
t3813 = -pkin(2) * t3875 + t4287;
t4108 = -0.2e1 * t4277;
t4119 = 0.2e1 * t3741 * t3720 * t4206 + pkin(2) * t3768;
t4137 = (t3901 + t3893) * t3660 + t3813 * t4298 + pkin(2) * t4349 + t4119 * qJ(3,2) + (-(-t3848 * t3875 + t3868) * t3881 + (t4108 + t4310) * pkin(2)) * t3768;
t3814 = -pkin(2) * t3877 + t4293;
t4107 = -0.2e1 * t4276;
t4118 = 0.2e1 * t3742 * t3721 * t4201 + pkin(2) * t3769;
t4136 = (t3901 + t3896) * t3661 + t3814 * t4297 + pkin(2) * t4350 + t4118 * qJ(3,1) + (-(-t3849 * t3877 + t3869) * t3883 + (t4107 + t4309) * pkin(2)) * t3769;
t4071 = t3737 * t4211;
t4095 = pkin(5) * t4270;
t4135 = t4095 - qJ(3,3) * (t3767 + t4071) + (t4278 - t4335) * t3767 - t4356;
t4069 = t3738 * t4206;
t4094 = pkin(5) * t4265;
t4134 = t4094 - qJ(3,2) * (t3768 + t4069) + (t4277 - t4336) * t3768 - t4355;
t4067 = t3739 * t4201;
t4093 = pkin(5) * t4260;
t4133 = t4093 - qJ(3,1) * (t3769 + t4067) + (t4276 - t4337) * t3769 - t4354;
t4014 = pkin(5) * t4071;
t4100 = pkin(5) * t4274;
t4132 = t4100 + (0.2e1 * t4335 + (-0.2e1 * t3862 + 0.2e1) * qJ(3,3)) * t3674 - t3873 * t4014 + (t4020 * t4199 + (t3701 * t4315 + t4321) * t3794) * t3773;
t4013 = pkin(5) * t4069;
t4098 = pkin(5) * t4273;
t4131 = t4098 + (0.2e1 * t4336 + (-0.2e1 * t3863 + 0.2e1) * qJ(3,2)) * t3675 - t3875 * t4013 + (t4019 * t4198 + (t3702 * t4314 + t4320) * t3797) * t3774;
t4012 = pkin(5) * t4067;
t4096 = pkin(5) * t4272;
t4130 = t4096 + (0.2e1 * t4337 + (-0.2e1 * t3864 + 0.2e1) * qJ(3,1)) * t3676 - t3877 * t4012 + (t4018 * t4197 + (t3703 * t4313 + t4319) * t3800) * t3775;
t4101 = pkin(5) * t3659 * t3873;
t4129 = -t4101 + 0.2e1 * (t4305 + t4340) * t3674 + (-0.4e1 * t3701 * t4052 - t4014) * t3879 + (-0.2e1 * t3826 * t3891 + t4315) * t3740 * t4228;
t4099 = pkin(5) * t3660 * t3875;
t4128 = -t4099 + 0.2e1 * (t4304 + t4339) * t3675 + (-0.4e1 * t3702 * t4049 - t4013) * t3881 + (-0.2e1 * t3827 * t3894 + t4314) * t3741 * t4225;
t4097 = pkin(5) * t3661 * t3877;
t4127 = -t4097 + 0.2e1 * (t4303 + t4338) * t3676 + (-0.4e1 * t3703 * t4046 - t4012) * t3883 + (-0.2e1 * t3828 * t3897 + t4313) * t3742 * t4222;
t4126 = pkin(5) * t4268 + 0.2e1 * qJ(3,3) * t3659 + (-t4199 - 0.2e1 * t4305) * t3767 + t4120;
t4125 = pkin(5) * t4263 + 0.2e1 * qJ(3,2) * t3660 + (-t4198 - 0.2e1 * t4304) * t3768 + t4119;
t4124 = pkin(5) * t4258 + 0.2e1 * qJ(3,1) * t3661 + (-t4197 - 0.2e1 * t4303) * t3769 + t4118;
t4112 = qJ(3,1) * t4312;
t4111 = qJ(3,2) * t4312;
t4110 = qJ(3,3) * t4312;
t4106 = -0.2e1 * t4270;
t4105 = -0.2e1 * t4265;
t4104 = -0.2e1 * t4260;
t4092 = pkin(5) * t4248;
t4091 = pkin(5) * t4245;
t4090 = pkin(5) * t4242;
t4083 = t3659 * t4213;
t4082 = t3660 * t4208;
t4081 = t3661 * t4203;
t4077 = t3794 * t4269;
t4076 = t3880 * t4271;
t4075 = t3797 * t4264;
t4074 = t3882 * t4266;
t4073 = t3800 * t4259;
t4072 = t3884 * t4261;
t4070 = t3796 * t4253;
t4068 = t3799 * t4252;
t4066 = t3802 * t4251;
t4065 = t3874 * t4250;
t4062 = t3876 * t4247;
t4059 = t3878 * t4244;
t4056 = t3770 * t4210;
t4055 = t3771 * t4205;
t4054 = t3772 * t4200;
t4044 = t3862 * t4219;
t4043 = t3863 * t4217;
t4042 = t3864 * t4215;
t4039 = t3796 * t4164;
t4036 = t3799 * t4160;
t4033 = t3802 * t4156;
t4008 = t3674 * t4044;
t3859 = t3873 ^ 2;
t4007 = t3859 * t4077;
t4006 = t4173 * t4271;
t4005 = t3674 * t4041;
t4004 = t4175 * t4269;
t4003 = t3675 * t4043;
t3860 = t3875 ^ 2;
t4002 = t3860 * t4075;
t4001 = t4170 * t4266;
t4000 = t3675 * t4038;
t3999 = t4172 * t4264;
t3998 = t3676 * t4042;
t3861 = t3877 ^ 2;
t3997 = t3861 * t4073;
t3996 = t4167 * t4261;
t3995 = t3676 * t4035;
t3994 = t4169 * t4259;
t3993 = (t3705 * t3879 + t3740 * t4214) * t4063;
t3992 = (t3707 * t3881 + t3741 * t4209) * t4060;
t3991 = (t3709 * t3883 + t3742 * t4204) * t4057;
t3990 = t3874 * t4070;
t3989 = t3880 * t4070;
t3988 = t3876 * t4068;
t3987 = t3882 * t4068;
t3986 = t3878 * t4066;
t3985 = t3884 * t4066;
t3981 = t3770 * t4040;
t3980 = t3770 * t4039;
t3979 = t3771 * t4037;
t3978 = t3771 * t4036;
t3977 = t3772 * t4034;
t3976 = t3772 * t4033;
t3975 = t3776 * t4041;
t3974 = t3777 * t4038;
t3973 = t3778 * t4035;
t3972 = t3873 * t4039;
t3971 = t3875 * t4036;
t3970 = t3877 * t4033;
t3965 = t3874 * t3993;
t3964 = t3876 * t3992;
t3963 = t3878 * t3991;
t3962 = t3873 * t3990;
t3961 = t3879 * t3990;
t3960 = t3875 * t3988;
t3959 = t3881 * t3988;
t3958 = t3877 * t3986;
t3957 = t3883 * t3986;
t3956 = t3835 * t3984;
t3955 = t3836 * t3983;
t3954 = t3837 * t3982;
t3935 = t3874 * t3956;
t3934 = t3876 * t3955;
t3933 = t3878 * t3954;
t3932 = -(t4268 * t4347 - t4101) * t3794 + (t3773 * t4311 + t3879 * t4092 / 0.2e1) * t4343;
t3931 = -(t4263 * t4347 - t4099) * t3797 + (t3774 * t4310 + t3881 * t4091 / 0.2e1) * t4342;
t3930 = -(t4258 * t4347 - t4097) * t3800 + (t3775 * t4309 + t3883 * t4090 / 0.2e1) * t4341;
t3929 = -(pkin(1) * t4106 - t4100) * t3794 + (t3879 * pkin(1) * t3773 - t3873 * t4092 / 0.2e1) * t4343;
t3928 = -(pkin(1) * t4105 - t4098) * t3797 + (t3881 * pkin(1) * t3774 - t3875 * t4091 / 0.2e1) * t4342;
t3927 = -(pkin(1) * t4104 - t4096) * t3800 + (t3883 * pkin(1) * t3775 - t3877 * t4090 / 0.2e1) * t4341;
t3926 = t3794 * (t3659 * t4196 + t3674 * t4241);
t3925 = t3794 * (-t3659 * t4187 + t3674 * t4240);
t3924 = t3797 * (t3660 * t4193 + t3675 * t4239);
t3923 = t3797 * (-t3660 * t4184 + t3675 * t4238);
t3922 = t3800 * (t3661 * t4190 + t3676 * t4237);
t3921 = t3800 * (-t3661 * t4181 + t3676 * t4236);
t3790 = (pkin(1) + 0.2e1 * t4300) * t3877 + t4107 + qJ(3,1);
t3789 = (pkin(1) + 0.2e1 * t4301) * t3875 + t4108 + qJ(3,2);
t3788 = (pkin(1) + 0.2e1 * t4302) * t3873 + t4109 + qJ(3,3);
t3786 = 0.2e1 * t4086 + t4115;
t3784 = 0.2e1 * t4085 + t4116;
t3782 = 0.2e1 * t4084 + t4117;
t3781 = qJ(3,1) * t3884 + t4115 * t3877;
t3780 = qJ(3,2) * t3882 + t4116 * t3875;
t3779 = qJ(3,3) * t3880 + t4117 * t3873;
t3751 = -0.2e1 * t3801 * t4230 + t3769;
t3750 = -0.2e1 * t3798 * t4231 + t3768;
t3749 = -0.2e1 * t3795 * t4232 + t3767;
t3748 = (t3843 * t4112 + t3846 * t4024) * t3864 + (t3786 * t4179 + t3843 * t4331) * t3883 + t3781 * t4294 + t3828 * t4188;
t3747 = (-t3843 * t4024 + t3846 * t4112) * t3864 + (-t3786 * t4188 + t3846 * t4331) * t3883 - t3781 * t4295 + t3828 * t4179;
t3746 = (t3842 * t4111 + t3845 * t4027) * t3863 + (t3784 * t4182 + t3842 * t4330) * t3881 + t3780 * t4288 + t3827 * t4191;
t3745 = (-t3842 * t4027 + t3845 * t4111) * t3863 + (-t3784 * t4191 + t3845 * t4330) * t3881 - t3780 * t4289 + t3827 * t4182;
t3744 = (t3841 * t4110 + t3844 * t4030) * t3862 + (t3782 * t4185 + t3841 * t4329) * t3879 + t3779 * t4282 + t3826 * t4194;
t3743 = (-t3841 * t4030 + t3844 * t4110) * t3862 + (-t3782 * t4194 + t3844 * t4329) * t3879 - t3779 * t4283 + t3826 * t4185;
t3724 = -t3769 + (t4230 - t4251) * t3801;
t3723 = -t3768 + (t4231 - t4252) * t3798;
t3722 = -t3767 + (t4232 - t4253) * t3795;
t3640 = t4350 - t4093 + 0.2e1 * t4306;
t3639 = t4349 - t4094 + 0.2e1 * t4307;
t3638 = t4348 - t4095 + 0.2e1 * t4308;
t3634 = t3938 + 0.2e1 * t4297;
t3633 = t3937 + 0.2e1 * t4298;
t3632 = t3936 + 0.2e1 * t4299;
t1 = [(-t3844 * t4077 - t3845 * t4075 - t3846 * t4073) * MDP(1) + (-t3844 * t4007 - t3845 * t4002 - t3846 * t3997 + (-t3757 * t3772 + t4059 * t4316) * t3970 + (-t3755 * t3771 + t4062 * t4317) * t3971 + (-t3753 * t3770 + t4065 * t4318) * t3972) * MDP(4) + (t3935 * t4318 + t3934 * t4317 + t3933 * t4316 + (t3751 * t4236 + t3994 * t4316) * t3800 + (t3750 * t4238 + t3999 * t4317) * t3797 + (t3749 * t4240 + t4004 * t4318) * t3794) * MDP(5) + (-t3844 * t3961 - t3845 * t3959 - t3846 * t3957 + t3873 * t3925 + t3875 * t3923 + t3877 * t3921) * MDP(6) + (t3844 * t3962 + t3845 * t3960 + t3846 * t3958 + t3879 * t3925 + t3881 * t3923 + t3883 * t3921) * MDP(7) + (t3753 * t4083 + t3755 * t4082 + t3757 * t4081) * MDP(8) + (t3930 * t4181 + t3931 * t4184 + t3932 * t4187 + (-t3753 * t4006 - t3755 * t4001 - t3757 * t3996) * pkin(5) + (t3753 * t3981 + t3755 * t3979 + t3757 * t3977) * pkin(1)) * MDP(9) + (t3927 * t4181 + t3928 * t4184 + t3929 * t4187 + (-t3753 * t4005 - t3755 * t4000 - t3757 * t3995) * pkin(5) + (t3753 * t3980 + t3755 * t3978 + t3757 * t3976) * pkin(1)) * MDP(10) + ((-t3748 * t4169 + t3757 * t3790) * t4054 + (-t3746 * t4172 + t3755 * t3789) * t4055 + (-t3744 * t4175 + t3753 * t3788) * t4056 + ((t3640 * t3757 - t3661 * t3748) * t3897 - t4127 * t4181) * t3800 + ((t3639 * t3755 - t3660 * t3746) * t3894 - t4128 * t4184) * t3797 + ((t3638 * t3753 - t3659 * t3744) * t3891 - t4129 * t4187) * t3794) * MDP(11) + (t3844 * t3965 + t3845 * t3964 + t3846 * t3963 + (-t3634 * t4181 + (t3748 * t3877 + t3757 * t3814) * t4257) * t3800 + (-t3633 * t4184 + (t3746 * t3875 + t3755 * t3813) * t4262) * t3797 + (-t3632 * t4187 + (t3744 * t3873 + t3753 * t3812) * t4267) * t3794) * MDP(12) + ((-t4130 * t4181 + (t3724 * t3748 + t4124 * t3757) * t3897) * t3800 + (-t4131 * t4184 + (t3723 * t3746 + t4125 * t3755) * t3894) * t3797 + (-t4132 * t4187 + (t3722 * t3744 + t4126 * t3753) * t3891) * t3794) * MDP(13) + ((-t4139 * t4181 + (t4133 * t3748 + t4136 * t3757) * t3897) * t3800 + (-t4140 * t4184 + (t4134 * t3746 + t4137 * t3755) * t3894) * t3797 + (-t4141 * t4187 + (t4135 * t3744 + t4138 * t3753) * t3891) * t3794) * MDP(14); (t3841 * t4077 + t3842 * t4075 + t3843 * t4073) * MDP(1) + (t3841 * t4007 + t3842 * t4002 + t3843 * t3997 + (-t3756 * t3772 + t4059 * t4344) * t3970 + (-t3754 * t3771 + t4062 * t4345) * t3971 + (-t3752 * t3770 + t4065 * t4346) * t3972) * MDP(4) + (t3935 * t4346 + t3934 * t4345 + t3933 * t4344 + (t3751 * t4237 + t3994 * t4344) * t3800 + (t3750 * t4239 + t3999 * t4345) * t3797 + (t3749 * t4241 + t4004 * t4346) * t3794) * MDP(5) + (t3841 * t3961 + t3842 * t3959 + t3843 * t3957 + t3873 * t3926 + t3875 * t3924 + t3877 * t3922) * MDP(6) + (-t3841 * t3962 - t3842 * t3960 - t3843 * t3958 + t3879 * t3926 + t3881 * t3924 + t3883 * t3922) * MDP(7) + (t3752 * t4083 + t3754 * t4082 + t3756 * t4081) * MDP(8) + (-t3930 * t4190 - t3931 * t4193 - t3932 * t4196 + (-t3752 * t4006 - t3754 * t4001 - t3756 * t3996) * pkin(5) + (t3752 * t3981 + t3754 * t3979 + t3756 * t3977) * pkin(1)) * MDP(9) + (-t3927 * t4190 - t3928 * t4193 - t3929 * t4196 + (-t3752 * t4005 - t3754 * t4000 - t3756 * t3995) * pkin(5) + (t3752 * t3980 + t3754 * t3978 + t3756 * t3976) * pkin(1)) * MDP(10) + ((-t3747 * t4169 + t3756 * t3790) * t4054 + (-t3745 * t4172 + t3754 * t3789) * t4055 + (-t3743 * t4175 + t3752 * t3788) * t4056 + ((t3640 * t3756 - t3661 * t3747) * t3897 + t4127 * t4190) * t3800 + ((t3639 * t3754 - t3660 * t3745) * t3894 + t4128 * t4193) * t3797 + ((t3638 * t3752 - t3659 * t3743) * t3891 + t4129 * t4196) * t3794) * MDP(11) + (-t3841 * t3965 - t3842 * t3964 - t3843 * t3963 + (t3634 * t4190 + (t3747 * t3877 + t3756 * t3814) * t4257) * t3800 + (t3633 * t4193 + (t3745 * t3875 + t3754 * t3813) * t4262) * t3797 + (t3632 * t4196 + (t3743 * t3873 + t3752 * t3812) * t4267) * t3794) * MDP(12) + ((t4130 * t4190 + (t3724 * t3747 + t3756 * t4124) * t3897) * t3800 + (t4131 * t4193 + (t3723 * t3745 + t3754 * t4125) * t3894) * t3797 + (t4132 * t4196 + (t3722 * t3743 + t3752 * t4126) * t3891) * t3794) * MDP(13) + ((t4139 * t4190 + (t3747 * t4133 + t3756 * t4136) * t3897) * t3800 + (t4140 * t4193 + (t3745 * t4134 + t3754 * t4137) * t3894) * t3797 + (t4141 * t4196 + (t3743 * t4135 + t3752 * t4138) * t3891) * t3794) * MDP(14); (-t4072 - t4074 - t4076) * MDP(1) + (-t3859 * t4076 - t3860 * t4074 - t3861 * t4072 + (t3778 * t4230 - 0.2e1 * t4158 * t4244) * t4034 + (t3777 * t4231 - 0.2e1 * t4162 * t4247) * t4037 + (t3776 * t4232 - 0.2e1 * t4166 * t4250) * t4040) * MDP(4) + (-0.2e1 * t3880 * t3956 - 0.2e1 * t3882 * t3955 - 0.2e1 * t3884 * t3954 + (-t3751 * t4215 + t3884 * t4104) * t3883 * t3800 + (-t3750 * t4217 + t3882 * t4105) * t3881 * t3797 + (-t3749 * t4219 + t3880 * t4106) * t3879 * t3794) * MDP(5) + (-t3879 * t3989 - t3881 * t3987 - t3883 * t3985 + (-t3676 * t3778 * t4156 - t3661 * t3884) * t4204 + (-t3675 * t3777 * t4160 - t3660 * t3882) * t4209 + (-t3674 * t3776 * t4164 - t3659 * t3880) * t4214) * MDP(6) + (t3873 * t3989 + t3875 * t3987 + t3877 * t3985 + (-t3661 * t4158 - t3998) * t3800 + (-t3660 * t4162 - t4003) * t3797 + (-t3659 * t4166 - t4008) * t3794) * MDP(7) + (-t3659 * t3975 - t3660 * t3974 - t3661 * t3973) * MDP(8) + (t3930 * t3884 + t3931 * t3882 + t3932 * t3880 + (t3973 * t4260 + t3974 * t4265 + t3975 * t4270) * pkin(5) + (-t3770 * t3776 * t3972 - t3771 * t3777 * t3971 - t3772 * t3778 * t3970) * pkin(1)) * MDP(9) + (t3927 * t3884 + t3928 * t3882 + t3929 * t3880 + (t3794 * t4008 + t3797 * t4003 + t3800 * t3998) * pkin(5) + (-t3770 * t3796 * t4044 - t3771 * t3799 * t4043 - t3772 * t3802 * t4042) * pkin(1)) * MDP(10) + ((-t3778 * t3790 - t4233) * t3976 + (-t3777 * t3789 - t4234) * t3978 + (-t3776 * t3788 - t4235) * t3980 + ((-t3640 * t4216 - t3661 * t3760) * t3897 - t4127 * t3884) * t3800 + ((-t3639 * t4218 - t3660 * t3759) * t3894 - t4128 * t3882) * t3797 + ((-t3638 * t4220 - t3659 * t3758) * t3891 - t4129 * t3880) * t3794) * MDP(11) + (t3880 * t3993 + t3882 * t3992 + t3884 * t3991 + (-t3884 * t3634 + (-t3814 * t4216 + t4233) * t4257) * t3800 + (-t3882 * t3633 + (-t3813 * t4218 + t4234) * t4262) * t3797 + (-t3880 * t3632 + (-t3812 * t4220 + t4235) * t4267) * t3794) * MDP(12) + ((-t4130 * t3884 + (t3724 * t3760 - t4124 * t4216) * t3897) * t3800 + (-t4131 * t3882 + (t3723 * t3759 - t4125 * t4218) * t3894) * t3797 + (-t4132 * t3880 + (t3722 * t3758 - t4126 * t4220) * t3891) * t3794) * MDP(13) + ((-t4139 * t3884 + (t3760 * t4133 - t4136 * t4216) * t3897) * t3800 + (-t4140 * t3882 + (t3759 * t4134 - t4137 * t4218) * t3894) * t3797 + (-t4141 * t3880 + (t3758 * t4135 - t4138 * t4220) * t3891) * t3794) * MDP(14);];
taucX  = t1;
