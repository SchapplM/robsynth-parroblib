% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRR1G3P3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR1G3P3A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR1G3P3A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:06:59
% EndTime: 2020-03-09 21:07:01
% DurationCPUTime: 2.50s
% Computational Cost: add. (22896->249), mult. (13391->428), div. (3141->10), fcn. (13530->36), ass. (0->209)
t979 = legFrame(1,2);
t970 = sin(t979);
t973 = cos(t979);
t990 = xDP(2);
t991 = xDP(1);
t1009 = t970 * t990 - t973 * t991;
t976 = pkin(7) + qJ(2,1);
t967 = qJ(3,1) + t976;
t955 = sin(t967);
t958 = cos(t967);
t989 = xDP(3);
t868 = t1009 * t958 + t955 * t989;
t978 = legFrame(2,2);
t969 = sin(t978);
t972 = cos(t978);
t1010 = t969 * t990 - t972 * t991;
t975 = pkin(7) + qJ(2,2);
t966 = qJ(3,2) + t975;
t954 = sin(t966);
t957 = cos(t966);
t867 = t1010 * t957 + t954 * t989;
t977 = legFrame(3,2);
t968 = sin(t977);
t971 = cos(t977);
t1011 = t968 * t990 - t971 * t991;
t974 = pkin(7) + qJ(2,3);
t965 = qJ(3,3) + t974;
t953 = sin(t965);
t956 = cos(t965);
t866 = t1011 * t956 + t953 * t989;
t995 = 1 / pkin(2);
t1097 = -2 * t995;
t961 = sin(t976);
t964 = cos(t976);
t889 = t955 * t964 - t961 * t958;
t1088 = 0.1e1 / t889;
t1096 = t1088 * t973;
t960 = sin(t975);
t963 = cos(t975);
t888 = t954 * t963 - t960 * t957;
t1089 = 0.1e1 / t888;
t1095 = t1089 * t972;
t959 = sin(t974);
t962 = cos(t974);
t887 = t953 * t962 - t959 * t956;
t1090 = 0.1e1 / t887;
t1094 = t1090 * t971;
t850 = -(-t1009 * t964 - t961 * t989) * pkin(2) + t868 * pkin(3);
t996 = 0.1e1 / pkin(2) ^ 2;
t1069 = t850 * t996;
t993 = 0.1e1 / pkin(3);
t1060 = t1088 * t993;
t1045 = t850 * t1060;
t1066 = t868 * t1088;
t847 = (t1045 - t1066) * t995;
t1021 = t847 * t1088 * t1069;
t1051 = t958 * t970;
t980 = xDDP(3);
t981 = xDDP(2);
t1093 = ((-t981 * t1051 - t955 * t980) * t995 + t1021) * t1088;
t849 = -(-t1010 * t963 - t960 * t989) * pkin(2) + t867 * pkin(3);
t1070 = t849 * t996;
t1062 = t1089 * t993;
t1046 = t849 * t1062;
t1067 = t867 * t1089;
t846 = (t1046 - t1067) * t995;
t1022 = t846 * t1089 * t1070;
t1052 = t957 * t969;
t1092 = ((-t981 * t1052 - t954 * t980) * t995 + t1022) * t1089;
t848 = -(-t1011 * t962 - t959 * t989) * pkin(2) + t866 * pkin(3);
t1071 = t848 * t996;
t1064 = t1090 * t993;
t1047 = t848 * t1064;
t1068 = t866 * t1090;
t845 = (t1047 - t1068) * t995;
t1023 = t845 * t1090 * t1071;
t1053 = t956 * t968;
t1091 = ((-t981 * t1053 - t953 * t980) * t995 + t1023) * t1090;
t882 = 0.1e1 / t887 ^ 2;
t884 = 0.1e1 / t888 ^ 2;
t886 = 0.1e1 / t889 ^ 2;
t1087 = g(1) / 0.2e1;
t1086 = -g(2) / 0.2e1;
t948 = -t977 + t965;
t1085 = sin(t948) / 0.2e1;
t950 = -t978 + t966;
t1084 = sin(t950) / 0.2e1;
t952 = -t979 + t967;
t1083 = sin(t952) / 0.2e1;
t947 = t977 + t965;
t1082 = cos(t947) / 0.2e1;
t949 = t978 + t966;
t1081 = cos(t949) / 0.2e1;
t951 = t979 + t967;
t1080 = cos(t951) / 0.2e1;
t1037 = t956 * t1094;
t982 = xDDP(1);
t1050 = t982 * t995;
t860 = t1037 * t1050;
t1002 = t860 + t1091;
t1041 = t866 * t882 * t996;
t1014 = t953 * t959 + t956 * t962;
t836 = t845 * pkin(3) - t1014 * t1068;
t830 = -t836 * t1041 + t1002;
t1079 = pkin(2) * t830;
t1035 = t957 * t1095;
t861 = t1035 * t1050;
t1001 = t861 + t1092;
t1040 = t867 * t884 * t996;
t1013 = t954 * t960 + t957 * t963;
t837 = t846 * pkin(3) - t1013 * t1067;
t831 = -t837 * t1040 + t1001;
t1078 = pkin(2) * t831;
t1033 = t958 * t1096;
t862 = t1033 * t1050;
t1000 = t862 + t1093;
t1039 = t868 * t886 * t996;
t1012 = t955 * t961 + t958 * t964;
t838 = t847 * pkin(3) - t1012 * t1066;
t832 = -t838 * t1039 + t1000;
t1077 = pkin(2) * t832;
t1076 = t981 - g(2);
t1075 = t982 - g(1);
t1048 = 0.2e1 * pkin(3);
t842 = (-t1068 + t1047 / 0.2e1) * t995;
t992 = pkin(3) ^ 2;
t1074 = (-t845 * t992 + (-t1014 * t842 * t1048 + t1068) * pkin(2)) * t993;
t843 = (-t1067 + t1046 / 0.2e1) * t995;
t1073 = (-t846 * t992 + (-t1013 * t843 * t1048 + t1067) * pkin(2)) * t993;
t844 = (-t1066 + t1045 / 0.2e1) * t995;
t1072 = (-t847 * t992 + (-t1012 * t844 * t1048 + t1066) * pkin(2)) * t993;
t1065 = t1090 * (pkin(2) * t959 + pkin(3) * t953);
t1063 = t1089 * (pkin(2) * t960 + pkin(3) * t954);
t1061 = t1088 * (pkin(2) * t961 + pkin(3) * t955);
t1059 = t1090 * t953;
t1058 = t1089 * t954;
t1057 = t1088 * t955;
t1049 = t993 * t995;
t1044 = t866 ^ 2 * t882 * t995;
t1043 = t867 ^ 2 * t884 * t995;
t1042 = t868 ^ 2 * t886 * t995;
t899 = pkin(2) * t962 + pkin(3) * t956;
t1038 = t1090 * t899 * t968;
t900 = pkin(2) * t963 + pkin(3) * t957;
t1036 = t1089 * t900 * t969;
t901 = pkin(2) * t964 + pkin(3) * t958;
t1034 = t1088 * t901 * t970;
t1032 = t899 * t1094;
t1031 = t1090 * t1053;
t1030 = t900 * t1095;
t1029 = t1089 * t1052;
t1028 = t901 * t1096;
t1027 = t1088 * t1051;
t1026 = t980 * t1049;
t1025 = t981 * t1049;
t1024 = t982 * t1049;
t1020 = -(t1014 * pkin(2) + pkin(3)) * t1023 * t1064 - t1024 * t1032 + t1025 * t1038 + t1026 * t1065;
t1019 = -(t1013 * pkin(2) + pkin(3)) * t1022 * t1062 - t1024 * t1030 + t1025 * t1036 + t1026 * t1063;
t1018 = -(t1012 * pkin(2) + pkin(3)) * t1021 * t1060 - t1024 * t1028 + t1025 * t1034 + t1026 * t1061;
t932 = sin(t947);
t939 = cos(t948);
t1017 = g(1) * t1085 + g(2) * t1082 + g(3) * t956 + t939 * t1086 + t932 * t1087;
t934 = sin(t949);
t941 = cos(t950);
t1016 = g(1) * t1084 + g(2) * t1081 + g(3) * t957 + t941 * t1086 + t934 * t1087;
t936 = sin(t951);
t943 = cos(t952);
t1015 = g(1) * t1083 + g(2) * t1080 + g(3) * t958 + t943 * t1086 + t936 * t1087;
t1008 = g(1) * t1082 + g(2) * t1085 - g(3) * t953 + t932 * t1086 + t939 * t1087;
t1007 = g(1) * t1081 + g(2) * t1084 - g(3) * t954 + t934 * t1086 + t941 * t1087;
t1006 = g(1) * t1080 + g(2) * t1083 - g(3) * t955 + t936 * t1086 + t943 * t1087;
t994 = 0.1e1 / pkin(3) ^ 2;
t988 = cos(qJ(3,1));
t987 = cos(qJ(3,2));
t986 = cos(qJ(3,3));
t985 = sin(qJ(3,1));
t984 = sin(qJ(3,2));
t983 = sin(qJ(3,3));
t904 = t973 * g(1) - t970 * g(2);
t903 = t972 * g(1) - t969 * g(2);
t902 = t971 * g(1) - t968 * g(2);
t877 = -g(3) * t961 + t904 * t964;
t876 = g(3) * t964 + t904 * t961;
t875 = -g(3) * t960 + t903 * t963;
t874 = g(3) * t963 + t903 * t960;
t873 = -g(3) * t959 + t902 * t962;
t872 = g(3) * t962 + t902 * t959;
t871 = t1075 * t970 + t1076 * t973;
t870 = t1075 * t969 + t1076 * t972;
t869 = t1075 * t968 + t1076 * t971;
t829 = t988 * t1042 - t985 * t1077 + t1006;
t828 = t987 * t1043 - t984 * t1078 + t1007;
t827 = t986 * t1044 - t983 * t1079 + t1008;
t826 = t985 * t1042 + t988 * t1077 + t1015;
t825 = t984 * t1043 + t987 * t1078 + t1016;
t824 = t983 * t1044 + t986 * t1079 + t1017;
t823 = (-t838 - t1072) * t1039 + t1000 + t1018;
t822 = 0.2e1 * t862 + (-0.2e1 * t838 - t1072) * t1039 + 0.2e1 * t1093 + t1018;
t821 = (-t837 - t1073) * t1040 + t1001 + t1019;
t820 = 0.2e1 * t861 + (-0.2e1 * t837 - t1073) * t1040 + 0.2e1 * t1092 + t1019;
t819 = (-t836 - t1074) * t1041 + t1002 + t1020;
t818 = 0.2e1 * t860 + (-0.2e1 * t836 - t1074) * t1041 + 0.2e1 * t1091 + t1020;
t817 = -(t985 * t822 + (t994 * t886 * t850 - 0.2e1 * t1060 * t1066) * t988 * t1069) * pkin(2) + t1006;
t816 = -(t984 * t820 + (t994 * t884 * t849 - 0.2e1 * t1062 * t1067) * t987 * t1070) * pkin(2) + t1007;
t815 = -(t983 * t818 + (t994 * t882 * t848 - 0.2e1 * t1064 * t1068) * t986 * t1071) * pkin(2) + t1008;
t814 = pkin(2) * (t1045 * t1097 * t844 * t985 + t822 * t988) + t1015;
t813 = pkin(2) * (t1046 * t1097 * t843 * t984 + t820 * t987) + t1016;
t812 = pkin(2) * (t1047 * t1097 * t842 * t983 + t818 * t986) + t1017;
t1 = [(t968 * t869 + t969 * t870 + t970 * t871) * MDP(1) + t1075 * MDP(8) + ((t832 * t1033 + t831 * t1035 + t830 * t1037) * MDP(2) + (t876 * t1033 + t874 * t1035 + t872 * t1037) * MDP(3) + (t877 * t1033 + t875 * t1035 + t873 * t1037) * MDP(4) + (t823 * t1033 + t821 * t1035 + t819 * t1037) * MDP(5) + (t814 * t1033 + t813 * t1035 + t812 * t1037) * MDP(6) + (t817 * t1033 + t816 * t1035 + t815 * t1037) * MDP(7) + ((-t823 * t1028 - t821 * t1030 - t819 * t1032) * MDP(5) + (-t826 * t1028 - t825 * t1030 - t824 * t1032) * MDP(6) + (-t829 * t1028 - t828 * t1030 - t827 * t1032) * MDP(7)) * t993) * t995; (t971 * t869 + t972 * t870 + t973 * t871) * MDP(1) + t1076 * MDP(8) + ((-t832 * t1027 - t831 * t1029 - t830 * t1031) * MDP(2) + (-t876 * t1027 - t874 * t1029 - t872 * t1031) * MDP(3) + (-t877 * t1027 - t875 * t1029 - t873 * t1031) * MDP(4) + (-t823 * t1027 - t821 * t1029 - t819 * t1031) * MDP(5) + (-t814 * t1027 - t813 * t1029 - t812 * t1031) * MDP(6) + (-t817 * t1027 - t816 * t1029 - t815 * t1031) * MDP(7) + ((t823 * t1034 + t821 * t1036 + t819 * t1038) * MDP(5) + (t826 * t1034 + t825 * t1036 + t824 * t1038) * MDP(6) + (t829 * t1034 + t828 * t1036 + t827 * t1038) * MDP(7)) * t993) * t995; (t980 - g(3)) * MDP(8) + ((-t832 * t1057 - t831 * t1058 - t830 * t1059) * MDP(2) + (-t876 * t1057 - t874 * t1058 - t872 * t1059) * MDP(3) + (-t877 * t1057 - t875 * t1058 - t873 * t1059) * MDP(4) + (-t823 * t1057 - t821 * t1058 - t819 * t1059) * MDP(5) + (-t814 * t1057 - t813 * t1058 - t812 * t1059) * MDP(6) + (-t817 * t1057 - t816 * t1058 - t815 * t1059) * MDP(7) + ((t823 * t1061 + t821 * t1063 + t819 * t1065) * MDP(5) + (t826 * t1061 + t825 * t1063 + t824 * t1065) * MDP(6) + (t829 * t1061 + t828 * t1063 + t827 * t1065) * MDP(7)) * t993) * t995;];
tauX  = t1;
