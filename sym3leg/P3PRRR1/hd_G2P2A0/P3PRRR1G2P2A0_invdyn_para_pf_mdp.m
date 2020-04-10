% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRR1G2P2A0
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
%   see P3PRRR1G2P2A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR1G2P2A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:30
% EndTime: 2020-03-09 21:18:33
% DurationCPUTime: 2.42s
% Computational Cost: add. (23274->258), mult. (13445->458), div. (3195->10), fcn. (13746->36), ass. (0->215)
t991 = pkin(7) + qJ(2,1);
t982 = qJ(3,1) + t991;
t970 = sin(t982);
t973 = cos(t982);
t976 = sin(t991);
t979 = cos(t991);
t901 = -t979 * t970 + t973 * t976;
t1101 = 0.1e1 / t901;
t1068 = t1101 * (pkin(2) * t979 + pkin(3) * t973);
t994 = legFrame(1,2);
t985 = sin(t994);
t1106 = t1101 * t985;
t990 = pkin(7) + qJ(2,2);
t981 = qJ(3,2) + t990;
t969 = sin(t981);
t972 = cos(t981);
t975 = sin(t990);
t978 = cos(t990);
t900 = -t978 * t969 + t972 * t975;
t1102 = 0.1e1 / t900;
t1070 = t1102 * (pkin(2) * t978 + pkin(3) * t972);
t993 = legFrame(2,2);
t984 = sin(t993);
t1105 = t1102 * t984;
t989 = pkin(7) + qJ(2,3);
t980 = qJ(3,3) + t989;
t968 = sin(t980);
t971 = cos(t980);
t974 = sin(t989);
t977 = cos(t989);
t899 = -t977 * t968 + t971 * t974;
t1103 = 0.1e1 / t899;
t1072 = t1103 * (pkin(2) * t977 + pkin(3) * t971);
t992 = legFrame(3,2);
t983 = sin(t992);
t1104 = t1103 * t983;
t1008 = 0.1e1 / pkin(3);
t1100 = 0.2e1 * t1008;
t1099 = -g(1) / 0.2e1;
t1098 = g(1) / 0.2e1;
t1097 = -g(2) / 0.2e1;
t1096 = g(2) / 0.2e1;
t962 = t992 + t980;
t1095 = sin(t962) / 0.2e1;
t964 = t993 + t981;
t1094 = sin(t964) / 0.2e1;
t966 = t994 + t982;
t1093 = sin(t966) / 0.2e1;
t963 = -t992 + t980;
t1092 = -cos(t963) / 0.2e1;
t965 = -t993 + t981;
t1091 = -cos(t965) / 0.2e1;
t967 = -t994 + t982;
t1090 = -cos(t967) / 0.2e1;
t1011 = 0.1e1 / pkin(2) ^ 2;
t1004 = xDP(3);
t1005 = xDP(2);
t1006 = xDP(1);
t986 = cos(t992);
t905 = -t1005 * t983 + t1006 * t986;
t1017 = -t1004 * t971 - t905 * t968;
t1044 = t1103 * t968 * t986;
t1010 = 0.1e1 / pkin(2);
t997 = xDDP(1);
t1058 = t1010 * t997;
t863 = t1044 * t1058;
t1050 = t968 * t1104;
t996 = xDDP(2);
t1059 = t1010 * t996;
t866 = t1050 * t1059;
t995 = xDDP(3);
t1060 = t1010 * t995;
t1071 = t1103 * t971;
t875 = t1060 * t1071;
t1039 = -t863 + t866 - t875;
t1078 = t1103 ^ 2;
t854 = -pkin(2) * (t1004 * t977 + t905 * t974) + t1017 * pkin(3);
t1057 = t854 * t1008;
t1063 = t1010 * t1103;
t851 = (-t1017 + t1057) * t1063;
t1084 = t851 * t854;
t1020 = t968 * t974 + t971 * t977;
t1081 = t1017 * t1103;
t845 = -pkin(3) * t851 + t1020 * t1081;
t839 = -(t1017 * t845 + t1084) * t1011 * t1078 + t1039;
t1089 = pkin(2) * t839;
t987 = cos(t993);
t906 = -t1005 * t984 + t1006 * t987;
t1016 = -t1004 * t972 - t906 * t969;
t1042 = t1102 * t969 * t987;
t864 = t1042 * t1058;
t1048 = t969 * t1105;
t867 = t1048 * t1059;
t1069 = t1102 * t972;
t876 = t1060 * t1069;
t1038 = -t864 + t867 - t876;
t1076 = t1102 ^ 2;
t855 = -pkin(2) * (t1004 * t978 + t906 * t975) + t1016 * pkin(3);
t1056 = t855 * t1008;
t1062 = t1010 * t1102;
t852 = (-t1016 + t1056) * t1062;
t1083 = t852 * t855;
t1019 = t969 * t975 + t972 * t978;
t1080 = t1016 * t1102;
t846 = -pkin(3) * t852 + t1019 * t1080;
t840 = -(t1016 * t846 + t1083) * t1011 * t1076 + t1038;
t1088 = pkin(2) * t840;
t988 = cos(t994);
t907 = -t1005 * t985 + t1006 * t988;
t1015 = -t1004 * t973 - t907 * t970;
t1040 = t1101 * t970 * t988;
t865 = t1040 * t1058;
t1046 = t970 * t1106;
t868 = t1046 * t1059;
t1067 = t1101 * t973;
t877 = t1060 * t1067;
t1037 = -t865 + t868 - t877;
t1074 = t1101 ^ 2;
t856 = -pkin(2) * (t1004 * t979 + t907 * t976) + t1015 * pkin(3);
t1055 = t856 * t1008;
t1061 = t1010 * t1101;
t853 = (-t1015 + t1055) * t1061;
t1082 = t853 * t856;
t1018 = t970 * t976 + t973 * t979;
t1079 = t1015 * t1101;
t847 = -pkin(3) * t853 + t1018 * t1079;
t841 = -(t1015 * t847 + t1082) * t1011 * t1074 + t1037;
t1087 = pkin(2) * t841;
t1086 = t996 - g(2);
t1085 = t997 - g(1);
t908 = pkin(2) * t974 + pkin(3) * t968;
t1077 = t1103 * t908;
t909 = pkin(2) * t975 + pkin(3) * t969;
t1075 = t1102 * t909;
t910 = pkin(2) * t976 + pkin(3) * t970;
t1073 = t1101 * t910;
t1007 = pkin(3) ^ 2;
t1054 = 0.2e1 * pkin(3);
t848 = (-t1017 + t1057 / 0.2e1) * t1063;
t1066 = t1008 * (t1007 * t851 + (t1020 * t848 * t1054 - t1081) * pkin(2));
t849 = (-t1016 + t1056 / 0.2e1) * t1062;
t1065 = t1008 * (t1007 * t852 + (t1019 * t849 * t1054 - t1080) * pkin(2));
t850 = (-t1015 + t1055 / 0.2e1) * t1061;
t1064 = t1008 * (t1007 * t853 + (t1018 * t850 * t1054 - t1079) * pkin(2));
t1053 = t1008 * t1010;
t1052 = 0.2e1 * t1053;
t1051 = t986 * t1077;
t1049 = t987 * t1075;
t1047 = t988 * t1073;
t1045 = t908 * t1104;
t1043 = t909 * t1105;
t1041 = t910 * t1106;
t894 = 0.1e1 / t899 ^ 2;
t1036 = t1010 * t1017 ^ 2 * t894;
t896 = 0.1e1 / t900 ^ 2;
t1035 = t1010 * t1016 ^ 2 * t896;
t898 = 0.1e1 / t901 ^ 2;
t1034 = t1010 * t1015 ^ 2 * t898;
t1033 = t995 * t1053;
t945 = sin(t963);
t950 = cos(t962);
t1029 = g(1) * t1092 + g(2) * t1095 + g(3) * t968 + t945 * t1097 + t950 * t1099;
t947 = sin(t965);
t952 = cos(t964);
t1028 = g(1) * t1091 + g(2) * t1094 + g(3) * t969 + t947 * t1097 + t952 * t1099;
t949 = sin(t967);
t954 = cos(t966);
t1027 = g(1) * t1090 + g(2) * t1093 + g(3) * t970 + t949 * t1097 + t954 * t1099;
t1026 = g(1) * t1095 + g(2) * t1092 + g(3) * t971 + t950 * t1096 + t945 * t1098;
t1025 = g(1) * t1094 + g(2) * t1091 + g(3) * t972 + t952 * t1096 + t947 * t1098;
t1024 = g(1) * t1093 + g(2) * t1090 + g(3) * t973 + t954 * t1096 + t949 * t1098;
t1023 = t851 * (t1020 * pkin(2) + pkin(3)) * t894 * t1057;
t1022 = t852 * (t1019 * pkin(2) + pkin(3)) * t896 * t1056;
t1021 = t853 * (t1018 * pkin(2) + pkin(3)) * t898 * t1055;
t1014 = t1033 * t1072 + (-t983 * t996 + t986 * t997) * t1053 * t1077;
t1013 = t1033 * t1070 + (-t984 * t996 + t987 * t997) * t1053 * t1075;
t1012 = t1033 * t1068 + (-t985 * t996 + t988 * t997) * t1053 * t1073;
t1009 = 0.1e1 / pkin(3) ^ 2;
t1003 = cos(qJ(3,1));
t1002 = cos(qJ(3,2));
t1001 = cos(qJ(3,3));
t1000 = sin(qJ(3,1));
t999 = sin(qJ(3,2));
t998 = sin(qJ(3,3));
t916 = -g(1) * t988 + g(2) * t985;
t915 = -g(1) * t987 + g(2) * t984;
t914 = -g(1) * t986 + g(2) * t983;
t889 = g(3) * t979 - t916 * t976;
t888 = g(3) * t978 - t915 * t975;
t887 = g(3) * t977 - t914 * t974;
t886 = g(3) * t976 + t916 * t979;
t885 = g(3) * t975 + t915 * t978;
t884 = g(3) * t974 + t914 * t977;
t880 = t1085 * t985 + t1086 * t988;
t879 = t1085 * t984 + t1086 * t987;
t878 = t1085 * t983 + t1086 * t986;
t838 = t1000 * t1034 + t1003 * t1087 + t1027;
t837 = t1002 * t1088 + t999 * t1035 + t1028;
t836 = t1001 * t1089 + t998 * t1036 + t1029;
t835 = -t1000 * t1087 + t1003 * t1034 + t1024;
t834 = t1002 * t1035 - t999 * t1088 + t1025;
t833 = t1001 * t1036 - t998 * t1089 + t1026;
t832 = (t1021 - (t1082 - (-t847 - t1064) * t1015) * t1074) * t1011 + t1012 + t1037;
t831 = -0.2e1 * t865 + 0.2e1 * t868 - 0.2e1 * t877 + (t1021 - (0.2e1 * t1082 - (-0.2e1 * t847 - t1064) * t1015) * t1074) * t1011 + t1012;
t830 = (t1022 - (t1083 - (-t846 - t1065) * t1016) * t1076) * t1011 + t1013 + t1038;
t829 = -0.2e1 * t864 + 0.2e1 * t867 - 0.2e1 * t876 + (t1022 - (0.2e1 * t1083 - (-0.2e1 * t846 - t1065) * t1016) * t1076) * t1011 + t1013;
t828 = (t1023 - (t1084 - (-t845 - t1066) * t1017) * t1078) * t1011 + t1014 + t1039;
t827 = -0.2e1 * t863 + 0.2e1 * t866 - 0.2e1 * t875 + (t1023 - (0.2e1 * t1084 - (-0.2e1 * t845 - t1066) * t1017) * t1078) * t1011 + t1014;
t826 = -pkin(2) * (t1000 * t831 + (t1009 * t856 - t1015 * t1100) * t898 * t856 * t1003 * t1011) + t1024;
t825 = -pkin(2) * (t999 * t829 + (t1009 * t855 - t1016 * t1100) * t896 * t855 * t1002 * t1011) + t1025;
t824 = -pkin(2) * (t998 * t827 + (t1009 * t854 - t1017 * t1100) * t894 * t854 * t1001 * t1011) + t1026;
t823 = pkin(2) * (-t1000 * t1052 * t1101 * t850 * t856 + t1003 * t831) + t1027;
t822 = pkin(2) * (-t1052 * t1102 * t849 * t855 * t999 + t1002 * t829) + t1028;
t821 = pkin(2) * (-t1052 * t1103 * t848 * t854 * t998 + t827 * t1001) + t1029;
t1 = [(t878 * t983 + t879 * t984 + t880 * t985) * MDP(1) + t1085 * MDP(8) + ((-t841 * t1040 - t840 * t1042 - t839 * t1044) * MDP(2) + (-t886 * t1040 - t885 * t1042 - t884 * t1044) * MDP(3) + (-t889 * t1040 - t888 * t1042 - t887 * t1044) * MDP(4) + (-t832 * t1040 - t830 * t1042 - t828 * t1044) * MDP(5) + (-t823 * t1040 - t822 * t1042 - t821 * t1044) * MDP(6) + (-t826 * t1040 - t825 * t1042 - t824 * t1044) * MDP(7) + ((t832 * t1047 + t830 * t1049 + t828 * t1051) * MDP(5) + (t838 * t1047 + t837 * t1049 + t836 * t1051) * MDP(6) + (t835 * t1047 + t834 * t1049 + t833 * t1051) * MDP(7)) * t1008) * t1010; (t878 * t986 + t879 * t987 + t880 * t988) * MDP(1) + t1086 * MDP(8) + ((t841 * t1046 + t840 * t1048 + t839 * t1050) * MDP(2) + (t886 * t1046 + t885 * t1048 + t884 * t1050) * MDP(3) + (t889 * t1046 + t888 * t1048 + t887 * t1050) * MDP(4) + (t832 * t1046 + t830 * t1048 + t828 * t1050) * MDP(5) + (t823 * t1046 + t822 * t1048 + t821 * t1050) * MDP(6) + (t826 * t1046 + t825 * t1048 + t824 * t1050) * MDP(7) + ((-t832 * t1041 - t830 * t1043 - t828 * t1045) * MDP(5) + (-t838 * t1041 - t837 * t1043 - t836 * t1045) * MDP(6) + (-t835 * t1041 - t834 * t1043 - t833 * t1045) * MDP(7)) * t1008) * t1010; (t995 - g(3)) * MDP(8) + ((-t841 * t1067 - t840 * t1069 - t839 * t1071) * MDP(2) + (-t886 * t1067 - t885 * t1069 - t884 * t1071) * MDP(3) + (-t889 * t1067 - t888 * t1069 - t887 * t1071) * MDP(4) + (-t832 * t1067 - t830 * t1069 - t828 * t1071) * MDP(5) + (-t823 * t1067 - t822 * t1069 - t821 * t1071) * MDP(6) + (-t826 * t1067 - t825 * t1069 - t824 * t1071) * MDP(7) + ((t832 * t1068 + t830 * t1070 + t828 * t1072) * MDP(5) + (t838 * t1068 + t837 * t1070 + t836 * t1072) * MDP(6) + (t835 * t1068 + t834 * t1070 + t833 * t1072) * MDP(7)) * t1008) * t1010;];
tauX  = t1;
