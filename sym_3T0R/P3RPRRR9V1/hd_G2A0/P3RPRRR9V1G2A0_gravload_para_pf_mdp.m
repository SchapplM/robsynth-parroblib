% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR9V1G2A0
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
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR9V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR9V1G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:53:14
% EndTime: 2020-08-06 18:53:17
% DurationCPUTime: 2.73s
% Computational Cost: add. (873->185), mult. (1101->324), div. (153->10), fcn. (1095->26), ass. (0->157)
t3908 = sin(pkin(7));
t3909 = cos(pkin(7));
t4031 = MDP(4) * t3909 - MDP(5) * t3908 + MDP(2);
t4014 = pkin(5) + pkin(6);
t3900 = qJ(2,1) + t4014;
t3888 = 0.1e1 / t3900;
t3924 = cos(qJ(1,1));
t3979 = t3888 * t3924;
t3918 = sin(qJ(1,1));
t3912 = legFrame(1,2);
t3891 = sin(t3912);
t3894 = cos(t3912);
t3940 = g(1) * t3894 - g(2) * t3891;
t4017 = g(3) * t3918 - t3940 * t3924;
t4027 = t4017 * t3979;
t3904 = pkin(7) + qJ(3,1);
t3884 = cos(t3904);
t3878 = 0.1e1 / t3884;
t3984 = t3878 * t3888;
t4026 = t4017 * t3984;
t3899 = qJ(2,2) + t4014;
t3887 = 0.1e1 / t3899;
t3922 = cos(qJ(1,2));
t3980 = t3887 * t3922;
t3916 = sin(qJ(1,2));
t3911 = legFrame(2,2);
t3890 = sin(t3911);
t3893 = cos(t3911);
t3941 = g(1) * t3893 - g(2) * t3890;
t4018 = g(3) * t3916 - t3941 * t3922;
t4025 = t4018 * t3980;
t3903 = pkin(7) + qJ(3,2);
t3883 = cos(t3903);
t3877 = 0.1e1 / t3883;
t3987 = t3877 * t3887;
t4024 = t4018 * t3987;
t3898 = qJ(2,3) + t4014;
t3886 = 0.1e1 / t3898;
t3920 = cos(qJ(1,3));
t3981 = t3886 * t3920;
t3914 = sin(qJ(1,3));
t3910 = legFrame(3,2);
t3889 = sin(t3910);
t3892 = cos(t3910);
t3942 = g(1) * t3892 - g(2) * t3889;
t4019 = g(3) * t3914 - t3942 * t3920;
t4023 = t4019 * t3981;
t3902 = pkin(7) + qJ(3,3);
t3882 = cos(t3902);
t3876 = 0.1e1 / t3882;
t3990 = t3876 * t3886;
t4022 = t4019 * t3990;
t4021 = MDP(3) - MDP(6);
t4016 = 2 * pkin(3);
t4015 = 0.2e1 * t3909 ^ 2;
t3919 = cos(qJ(3,3));
t3905 = t3919 ^ 2;
t4010 = t3905 * pkin(3);
t3921 = cos(qJ(3,2));
t3906 = t3921 ^ 2;
t4009 = t3906 * pkin(3);
t3923 = cos(qJ(3,1));
t3907 = t3923 ^ 2;
t4008 = t3907 * pkin(3);
t4007 = t3919 * pkin(2);
t4006 = t3921 * pkin(2);
t4005 = t3923 * pkin(2);
t3895 = t3918 * pkin(1);
t3827 = -g(3) * (qJ(2,1) * t3924 - t3895) - t3940 * (pkin(1) * t3924 + qJ(2,1) * t3918);
t4002 = t3827 * t3878;
t3896 = pkin(1) * t3914;
t3828 = -g(3) * (qJ(2,3) * t3920 - t3896) - t3942 * (pkin(1) * t3920 + qJ(2,3) * t3914);
t4001 = t3828 * t3876;
t3897 = pkin(1) * t3916;
t3829 = -g(3) * (qJ(2,2) * t3922 - t3897) - t3941 * (pkin(1) * t3922 + qJ(2,2) * t3916);
t4000 = t3829 * t3877;
t3999 = t4019 * t3886;
t3998 = t4018 * t3887;
t3997 = t4017 * t3888;
t3913 = sin(qJ(3,3));
t3972 = t3908 * t3913;
t3996 = t4019 / (t3909 * t3919 - t3972);
t3915 = sin(qJ(3,2));
t3971 = t3908 * t3915;
t3995 = t4018 / (t3909 * t3921 - t3971);
t3917 = sin(qJ(3,1));
t3970 = t3908 * t3917;
t3994 = t4017 / (t3909 * t3923 - t3970);
t3926 = pkin(2) / 0.2e1;
t3993 = (pkin(3) * t3919 + t3926) * t3913;
t3992 = (pkin(3) * t3921 + t3926) * t3915;
t3991 = (pkin(3) * t3923 + t3926) * t3917;
t3989 = t3876 * t3889;
t3988 = t3876 * t3892;
t3986 = t3877 * t3890;
t3985 = t3877 * t3893;
t3983 = t3878 * t3891;
t3982 = t3878 * t3894;
t3978 = t3889 * t3914;
t3977 = t3890 * t3916;
t3976 = t3891 * t3918;
t3975 = t3892 * t3914;
t3974 = t3893 * t3916;
t3973 = t3894 * t3918;
t3885 = pkin(1) * t3908;
t3969 = t3919 * (-pkin(3) * t3913 + t3885);
t3968 = t3921 * (-pkin(3) * t3915 + t3885);
t3967 = t3923 * (-pkin(3) * t3917 + t3885);
t3936 = g(3) * t3920 + t3942 * t3914;
t3954 = t3936 * t3990;
t3935 = g(3) * t3922 + t3941 * t3916;
t3953 = t3935 * t3987;
t3934 = g(3) * t3924 + t3940 * t3918;
t3952 = t3934 * t3984;
t3951 = t3914 * t3972;
t3950 = t3916 * t3971;
t3949 = t3918 * t3970;
t3948 = -t3900 * t3924 + t3895;
t3947 = -t3898 * t3920 + t3896;
t3946 = -t3899 * t3922 + t3897;
t3879 = sin(t3902);
t3945 = t3879 * t4022;
t3880 = sin(t3903);
t3944 = t3880 * t4024;
t3881 = sin(t3904);
t3943 = t3881 * t4026;
t3939 = pkin(2) * t3951 + (t3951 * t4016 - t3947) * t3919;
t3938 = pkin(2) * t3950 + (t3950 * t4016 - t3946) * t3921;
t3937 = pkin(2) * t3949 + (t3949 * t4016 - t3948) * t3923;
t3927 = 1 / pkin(3);
t3925 = -pkin(3) / 0.2e1;
t3875 = pkin(2) * t3909 + pkin(1);
t3868 = t4008 + t4005 / 0.2e1 + t3925;
t3867 = t4009 + t4006 / 0.2e1 + t3925;
t3866 = t4010 + t4007 / 0.2e1 + t3925;
t3865 = g(1) * t3891 + g(2) * t3894;
t3864 = g(1) * t3890 + g(2) * t3893;
t3863 = g(1) * t3889 + g(2) * t3892;
t3841 = t3881 * t3891 + t3884 * t3973;
t3840 = t3881 * t3894 - t3884 * t3976;
t3839 = t3880 * t3890 + t3883 * t3974;
t3838 = t3880 * t3893 - t3883 * t3977;
t3837 = t3879 * t3889 + t3882 * t3975;
t3836 = t3879 * t3892 - t3882 * t3978;
t3835 = pkin(1) * t3917 + (-pkin(3) + t4005 + 0.2e1 * t4008) * t3908;
t3834 = pkin(1) * t3915 + (-pkin(3) + t4006 + 0.2e1 * t4009) * t3908;
t3833 = pkin(1) * t3913 + (-pkin(3) + t4007 + 0.2e1 * t4010) * t3908;
t3832 = t3948 * t3970 + (t3907 - 0.1e1) * t3918 * pkin(3);
t3831 = t3946 * t3971 + (t3906 - 0.1e1) * t3916 * pkin(3);
t3830 = t3947 * t3972 + (t3905 - 0.1e1) * t3914 * pkin(3);
t3826 = t3863 * t3879 + t3882 * t3936;
t3825 = -t3863 * t3882 + t3879 * t3936;
t3824 = t3865 * t3881 + t3884 * t3934;
t3823 = -t3865 * t3884 + t3881 * t3934;
t3822 = t3864 * t3880 + t3883 * t3935;
t3821 = -t3864 * t3883 + t3880 * t3935;
t1 = [((t3841 * t4002 - ((t3868 * t3973 + t3891 * t3991) * t4015 + (t3891 * t3835 - t3937 * t3894) * t3909 - t3832 * t3894 + t3891 * t3967) * t3994) * t3888 + (t3839 * t4000 - ((t3867 * t3974 + t3890 * t3992) * t4015 + (t3890 * t3834 - t3938 * t3893) * t3909 - t3831 * t3893 + t3890 * t3968) * t3995) * t3887 + (t3837 * t4001 - ((t3866 * t3975 + t3889 * t3993) * t4015 + (t3889 * t3833 - t3939 * t3892) * t3909 - t3830 * t3892 + t3889 * t3969) * t3996) * t3886) * MDP(7) + (t3837 * t3999 + t3839 * t3998 + t3841 * t3997 + (t3821 * t3986 + t3823 * t3983 + t3825 * t3989) * t3927) * MDP(13) + (-t3837 * t3945 - t3839 * t3944 - t3841 * t3943 + (t3822 * t3986 + t3824 * t3983 + t3826 * t3989) * t3927) * MDP(14) - g(1) * MDP(15) + t4021 * (t3837 * t3954 + t3839 * t3953 + t3841 * t3952) + t4031 * (t3837 * t4022 + t3839 * t4024 + t3841 * t4026); ((t3840 * t4002 - ((-t3868 * t3976 + t3894 * t3991) * t4015 + (t3894 * t3835 + t3937 * t3891) * t3909 + t3832 * t3891 + t3894 * t3967) * t3994) * t3888 + (t3838 * t4000 - ((-t3867 * t3977 + t3893 * t3992) * t4015 + (t3893 * t3834 + t3938 * t3890) * t3909 + t3831 * t3890 + t3893 * t3968) * t3995) * t3887 + (t3836 * t4001 - ((-t3866 * t3978 + t3892 * t3993) * t4015 + (t3892 * t3833 + t3939 * t3889) * t3909 + t3830 * t3889 + t3892 * t3969) * t3996) * t3886) * MDP(7) + (t3836 * t3999 + t3838 * t3998 + t3840 * t3997 + (t3821 * t3985 + t3823 * t3982 + t3825 * t3988) * t3927) * MDP(13) + (-t3836 * t3945 - t3838 * t3944 - t3840 * t3943 + (t3822 * t3985 + t3824 * t3982 + t3826 * t3988) * t3927) * MDP(14) - g(2) * MDP(15) + t4021 * (t3836 * t3954 + t3838 * t3953 + t3840 * t3952) + t4031 * (t3836 * t4022 + t3838 * t4024 + t3840 * t4026); ((-t3918 * t3900 * t4017 + (t3827 - (pkin(3) * t3884 + t3875) * t4017) * t3924) * t3888 + (-t3916 * t3899 * t4018 + (t3829 - (pkin(3) * t3883 + t3875) * t4018) * t3922) * t3887 + (-t3914 * t3898 * t4019 + (t3828 - (pkin(3) * t3882 + t3875) * t4019) * t3920) * t3886) * MDP(7) + (t3882 * t4023 + t3883 * t4025 + t3884 * t4027) * MDP(13) + (-t3879 * t4023 - t3880 * t4025 - t3881 * t4027) * MDP(14) - g(3) * MDP(15) + t4021 * (t3934 * t3979 + t3935 * t3980 + t3936 * t3981) + t4031 * (t4027 + t4025 + t4023);];
taugX  = t1;
