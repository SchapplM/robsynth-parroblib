% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRRRR1V1G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:33:33
% EndTime: 2020-08-07 03:33:40
% DurationCPUTime: 7.76s
% Computational Cost: add. (37809->459), mult. (68160->754), div. (16335->8), fcn. (72651->54), ass. (0->355)
t963 = qJ(2,3) + qJ(3,3);
t1190 = cos(qJ(1,3) - t963) + cos(qJ(1,3) + t963);
t965 = qJ(2,2) + qJ(3,2);
t1189 = cos(qJ(1,2) - t965) + cos(qJ(1,2) + t965);
t967 = qJ(2,1) + qJ(3,1);
t1188 = cos(qJ(1,1) - t967) + cos(qJ(1,1) + t967);
t1006 = (pkin(2) ^ 2);
t1000 = (rSges(3,2) ^ 2);
t1002 = (rSges(3,1) ^ 2);
t938 = t1000 + t1002;
t907 = (t1006 + t938);
t1001 = (rSges(2,2) ^ 2);
t1003 = (rSges(2,1) ^ 2);
t940 = (t1001 + t1003);
t996 = 2 * pkin(1) ^ 2;
t1008 = -((2 * rSges(3,3) ^ 2 + t907 + t996) * m(3)) / 0.2e1 - ((2 * rSges(2,3) ^ 2 + t940 + t996) * m(2)) / 0.2e1 - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - Icges(3,2) / 0.2e1 - Icges(2,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,1) / 0.2e1 - Icges(1,3);
t1171 = m(2) * rSges(2,2);
t1172 = pkin(2) * m(3);
t908 = m(2) * rSges(2,1) + t1172;
t978 = sin(qJ(2,1));
t987 = cos(qJ(2,1));
t1015 = t1171 * t987 + t908 * t978;
t1184 = 2 * pkin(1);
t1045 = m(3) * t1184;
t1046 = t1171 * t1184;
t925 = sin(t967);
t1160 = pkin(1) * t925;
t1071 = 0.2e1 * t1160;
t1074 = -2 * pkin(1) * t908;
t919 = rSges(2,1) * t1171 - Icges(2,4);
t999 = 0.2e1 * qJ(2,1);
t949 = cos(t999);
t1109 = t919 * t949;
t914 = 0.2e1 * t967;
t896 = cos(t914);
t1169 = m(3) * rSges(3,2);
t917 = rSges(3,1) * t1169 - Icges(3,4);
t1112 = t917 * t896;
t885 = (m(3) * (-t1000 + t1002)) - Icges(3,1) + Icges(3,2);
t893 = sin(t914);
t1115 = t885 * t893;
t871 = (t1006 * m(3)) + ((-t1001 + t1003) * m(2)) + Icges(2,2) - Icges(2,1);
t946 = sin(t999);
t1118 = t871 * t946;
t1005 = 1 / pkin(3);
t1007 = 1 / pkin(2);
t1077 = t1005 * t1007;
t993 = xDP(3);
t1040 = t993 * t1077;
t977 = sin(qJ(3,1));
t1094 = t977 * t978;
t1054 = pkin(3) * t1094;
t986 = cos(qJ(3,1));
t1139 = pkin(3) * t986;
t911 = pkin(2) + t1139;
t884 = t911 * t987;
t863 = t884 - t1054;
t952 = 0.1e1 / t977;
t988 = cos(qJ(1,1));
t1041 = t863 * t952 * t988;
t1012 = t1040 * t1041;
t995 = xDP(1);
t1038 = t995 * t1077;
t970 = legFrame(1,2);
t937 = cos(t970);
t979 = sin(qJ(1,1));
t1103 = t937 * t979;
t1093 = t977 * t987;
t1053 = pkin(3) * t1093;
t860 = t911 * t978 + t1053;
t934 = sin(t970);
t842 = -t1103 * t863 + t860 * t934;
t821 = t842 * t952 * t1038;
t994 = xDP(2);
t1039 = t994 * t1077;
t1106 = t934 * t979;
t845 = t1106 * t863 + t860 * t937;
t824 = t845 * t952 * t1039;
t809 = t821 + t824 - t1012;
t879 = rSges(3,1) * t977 + rSges(3,2) * t986;
t1127 = t809 * t879;
t966 = t999 + qJ(3,1);
t924 = sin(t966);
t1148 = pkin(2) * t924;
t931 = cos(t967);
t1157 = pkin(1) * t931;
t1166 = -t885 / 0.2e1;
t1167 = -t871 / 0.2e1;
t930 = cos(t966);
t1178 = rSges(3,1) * t924 + rSges(3,2) * t930;
t959 = t986 ^ 2;
t1142 = pkin(3) * t959;
t1009 = pkin(1) * t1094 - pkin(3) + t1142;
t1034 = (t1007 * t993) / 0.2e1;
t1080 = t1007 * t995;
t1081 = t1007 * t994;
t1121 = t1188 * t952;
t1018 = -t986 * t987 + t1094;
t1019 = -t978 * t986 - t1093;
t838 = t1018 * t1106 + t1019 * t937;
t839 = -t1018 * t1103 + t1019 * t934;
t812 = t1034 * t1121 + (t1080 * t839 + t1081 * t838) * t952;
t1154 = pkin(2) * t812;
t1047 = t1154 / 0.2e1;
t943 = t986 * pkin(2);
t1050 = t812 * t943;
t806 = t809 + t812;
t1059 = t806 * t1139;
t1082 = t1007 * t952;
t1091 = t979 * t988;
t1092 = t978 * t987;
t1130 = t806 * t809;
t1145 = pkin(3) * t806;
t1181 = 0.2e1 * t988 ^ 2;
t800 = t1059 + t1154;
t866 = 0.1e1 / (pkin(2) * t987 + pkin(3) * t931 + pkin(1));
t818 = (-t979 * t993 + (-t934 * t994 + t937 * t995) * t988) * t866;
t876 = -pkin(3) + t943 + 0.2e1 * t1142;
t889 = pkin(1) - 0.2e1 * t1054;
t960 = t987 ^ 2;
t776 = (-pkin(3) * t1130 + (-0.4e1 * ((t986 * t1047 + (t959 - 0.1e1 / 0.2e1) * t1145) * t1092 + ((t1047 + t1059) * t960 - t800 / 0.2e1) * t977) * t1091 + t1188 * ((t1053 * t806 + t800 * t978) * t979 - t988 * (pkin(1) + t863) * t818) + (t1181 - 0.2e1) * t818 * (t876 * t960 + (-pkin(2) * t1094 + t889 * t986) * t987 - t1009)) * t818 / 0.2e1 + (((t1050 + (0.2e1 * t959 - 0.1e1) * t1145) * t960 - (0.2e1 * t1059 + t1154) * t977 * t1092 + t1145 - t1145 * t959) * t1181 - 0.2e1 * t818 * (t876 * t1092 + ((pkin(2) + 0.2e1 * t1139) * t960 - t911) * t977) * t1091 - 0.2e1 * t1050 - 0.2e1 * t1145 + t1188 * ((t1054 * t806 - t800 * t987) * t988 + t979 * t860 * t818)) * t812 / 0.2e1) * t1082;
t1004 = pkin(3) ^ 2;
t1065 = pkin(3) * t943;
t1070 = t1004 - t1006;
t1174 = -0.2e1 * t1004;
t803 = t821 / 0.2e1 + t824 / 0.2e1 - t1012 / 0.2e1 + t812;
t815 = t818 ^ 2;
t785 = ((t943 + pkin(3)) * t1130 + (-t815 * ((t1174 * t959 - 0.2e1 * t1065 + t1070) * t960 - t889 * t884 + pkin(3) * t1009) + (t1004 * t806 + t1006 * t812 + 0.2e1 * t1065 * t803) * t812) * t1005) * t1082;
t788 = 0.2e1 * (-t1145 * t925 - t1154 * t978) * t866 * t818;
t1168 = m(3) * rSges(3,3);
t915 = rSges(3,2) * t1168 - Icges(3,6);
t916 = rSges(3,1) * t1168 - Icges(3,5);
t851 = -t915 * t931 - t916 * t925;
t1170 = m(2) * rSges(2,3);
t1069 = rSges(2,1) * t1170;
t1075 = pkin(2) * t1168;
t880 = -Icges(2,5) + t1069 + t1075;
t1068 = rSges(2,2) * t1170;
t918 = -Icges(2,6) + t1068;
t848 = -t880 * t978 - t918 * t987 + t851;
t872 = -t1069 / 0.2e1 - t1075 / 0.2e1 + Icges(2,5) / 0.2e1;
t897 = t1068 / 0.2e1 - Icges(2,6) / 0.2e1;
t1076 = pkin(2) * t1169;
t900 = t977 * t1076;
t1031 = t866 * (t848 * t776 + t851 * t785 + (t978 * t1046 + t987 * t1074 + t896 * t1166 + t949 * t1167 + t1008 + t900 + t919 * t946 + t917 * t893 + ((t1071 + t1148) * rSges(3,2) + (-0.2e1 * t1157 + (-t930 - t986) * pkin(2)) * rSges(3,1)) * m(3)) * t788 + (0.2e1 * (-t872 * t987 - t897 * t978) * t812 + (-t1015 * t1184 - 0.2e1 * t1109 - t1118) * t818) * t812 + (-0.2e1 * t1178 * t803 - t1127) * t818 * t1172 + ((-t1115 - 0.2e1 * t1112 + (-rSges(3,1) * t925 - rSges(3,2) * t931) * t1045) * t818 + (-t915 * t925 + t916 * t931) * t806) * t806);
t1187 = t1031 * t988;
t975 = sin(qJ(2,2));
t984 = cos(qJ(2,2));
t1016 = t1171 * t984 + t908 * t975;
t923 = sin(t965);
t1161 = pkin(1) * t923;
t1072 = 0.2e1 * t1161;
t998 = 0.2e1 * qJ(2,2);
t948 = cos(t998);
t1110 = t919 * t948;
t913 = 0.2e1 * t965;
t895 = cos(t913);
t1113 = t917 * t895;
t892 = sin(t913);
t1116 = t885 * t892;
t945 = sin(t998);
t1119 = t871 * t945;
t974 = sin(qJ(3,2));
t1098 = t974 * t975;
t1056 = pkin(3) * t1098;
t983 = cos(qJ(3,2));
t1140 = pkin(3) * t983;
t910 = pkin(2) + t1140;
t883 = t910 * t984;
t862 = t883 - t1056;
t951 = 0.1e1 / t974;
t985 = cos(qJ(1,2));
t1042 = t862 * t951 * t985;
t1013 = t1040 * t1042;
t969 = legFrame(2,2);
t936 = cos(t969);
t976 = sin(qJ(1,2));
t1104 = t936 * t976;
t1097 = t974 * t984;
t1055 = pkin(3) * t1097;
t859 = t910 * t975 + t1055;
t933 = sin(t969);
t841 = -t1104 * t862 + t859 * t933;
t820 = t841 * t951 * t1038;
t1107 = t933 * t976;
t844 = t1107 * t862 + t859 * t936;
t823 = t844 * t951 * t1039;
t808 = t820 + t823 - t1013;
t878 = rSges(3,1) * t974 + rSges(3,2) * t983;
t1128 = t808 * t878;
t964 = qJ(3,2) + t998;
t922 = sin(t964);
t1149 = pkin(2) * t922;
t929 = cos(t965);
t1158 = pkin(1) * t929;
t928 = cos(t964);
t1177 = rSges(3,1) * t922 + rSges(3,2) * t928;
t956 = t983 ^ 2;
t1143 = pkin(3) * t956;
t1010 = pkin(1) * t1098 - pkin(3) + t1143;
t1122 = t1189 * t951;
t1020 = -t983 * t984 + t1098;
t1021 = -t975 * t983 - t1097;
t836 = t1020 * t1107 + t1021 * t936;
t837 = -t1020 * t1104 + t1021 * t933;
t811 = t1034 * t1122 + (t1080 * t837 + t1081 * t836) * t951;
t1155 = pkin(2) * t811;
t1048 = t1155 / 0.2e1;
t942 = t983 * pkin(2);
t1051 = t811 * t942;
t805 = t808 + t811;
t1060 = t805 * t1140;
t1083 = t1007 * t951;
t1095 = t976 * t985;
t1096 = t975 * t984;
t1131 = t805 * t808;
t1146 = pkin(3) * t805;
t1182 = 0.2e1 * t985 ^ 2;
t799 = t1060 + t1155;
t865 = 0.1e1 / (pkin(2) * t984 + pkin(3) * t929 + pkin(1));
t817 = (-t976 * t993 + (-t933 * t994 + t936 * t995) * t985) * t865;
t875 = -pkin(3) + t942 + 0.2e1 * t1143;
t888 = pkin(1) - 0.2e1 * t1056;
t957 = t984 ^ 2;
t775 = (-pkin(3) * t1131 + (-0.4e1 * ((t983 * t1048 + (t956 - 0.1e1 / 0.2e1) * t1146) * t1096 + ((t1048 + t1060) * t957 - t799 / 0.2e1) * t974) * t1095 + t1189 * ((t1055 * t805 + t799 * t975) * t976 - t985 * (pkin(1) + t862) * t817) + (t1182 - 0.2e1) * t817 * (t875 * t957 + (-pkin(2) * t1098 + t888 * t983) * t984 - t1010)) * t817 / 0.2e1 + (((t1051 + (0.2e1 * t956 - 0.1e1) * t1146) * t957 - (0.2e1 * t1060 + t1155) * t974 * t1096 + t1146 - t1146 * t956) * t1182 - 0.2e1 * t817 * (t875 * t1096 + ((pkin(2) + 0.2e1 * t1140) * t957 - t910) * t974) * t1095 - 0.2e1 * t1051 - 0.2e1 * t1146 + t1189 * ((t1056 * t805 - t799 * t984) * t985 + t976 * t859 * t817)) * t811 / 0.2e1) * t1083;
t1066 = pkin(3) * t942;
t802 = t820 / 0.2e1 + t823 / 0.2e1 - t1013 / 0.2e1 + t811;
t814 = t817 ^ 2;
t784 = ((t942 + pkin(3)) * t1131 + (-t814 * ((t1174 * t956 - 0.2e1 * t1066 + t1070) * t957 - t888 * t883 + pkin(3) * t1010) + (t1004 * t805 + t1006 * t811 + 0.2e1 * t1066 * t802) * t811) * t1005) * t1083;
t787 = 0.2e1 * (-t1146 * t923 - t1155 * t975) * t865 * t817;
t850 = -t915 * t929 - t916 * t923;
t847 = -t880 * t975 - t918 * t984 + t850;
t899 = t974 * t1076;
t1032 = t865 * (t847 * t775 + t850 * t784 + (t975 * t1046 + t984 * t1074 + t895 * t1166 + t948 * t1167 + t1008 + t899 + t919 * t945 + t917 * t892 + ((t1072 + t1149) * rSges(3,2) + (-0.2e1 * t1158 + (-t928 - t983) * pkin(2)) * rSges(3,1)) * m(3)) * t787 + (0.2e1 * (-t872 * t984 - t897 * t975) * t811 + (-t1016 * t1184 - 0.2e1 * t1110 - t1119) * t817) * t811 + (-0.2e1 * t1177 * t802 - t1128) * t817 * t1172 + ((-t1116 - 0.2e1 * t1113 + (-rSges(3,1) * t923 - rSges(3,2) * t929) * t1045) * t817 + (-t915 * t923 + t916 * t929) * t805) * t805);
t1186 = t1032 * t985;
t972 = sin(qJ(2,3));
t981 = cos(qJ(2,3));
t1017 = t1171 * t981 + t908 * t972;
t921 = sin(t963);
t1162 = pkin(1) * t921;
t1073 = 0.2e1 * t1162;
t997 = 0.2e1 * qJ(2,3);
t947 = cos(t997);
t1111 = t919 * t947;
t912 = 0.2e1 * t963;
t894 = cos(t912);
t1114 = t917 * t894;
t891 = sin(t912);
t1117 = t885 * t891;
t944 = sin(t997);
t1120 = t871 * t944;
t971 = sin(qJ(3,3));
t1102 = t971 * t972;
t1058 = pkin(3) * t1102;
t980 = cos(qJ(3,3));
t1141 = pkin(3) * t980;
t909 = pkin(2) + t1141;
t882 = t909 * t981;
t861 = t882 - t1058;
t950 = 0.1e1 / t971;
t982 = cos(qJ(1,3));
t1043 = t861 * t950 * t982;
t1014 = t1040 * t1043;
t968 = legFrame(3,2);
t935 = cos(t968);
t973 = sin(qJ(1,3));
t1105 = t935 * t973;
t1101 = t971 * t981;
t1057 = pkin(3) * t1101;
t858 = t909 * t972 + t1057;
t932 = sin(t968);
t840 = -t1105 * t861 + t858 * t932;
t819 = t840 * t950 * t1038;
t1108 = t932 * t973;
t843 = t1108 * t861 + t858 * t935;
t822 = t843 * t950 * t1039;
t807 = t819 + t822 - t1014;
t877 = rSges(3,1) * t971 + rSges(3,2) * t980;
t1129 = t807 * t877;
t962 = t997 + qJ(3,3);
t920 = sin(t962);
t1150 = pkin(2) * t920;
t927 = cos(t963);
t1159 = pkin(1) * t927;
t926 = cos(t962);
t1176 = rSges(3,1) * t920 + rSges(3,2) * t926;
t953 = t980 ^ 2;
t1144 = pkin(3) * t953;
t1011 = pkin(1) * t1102 - pkin(3) + t1144;
t1123 = t1190 * t950;
t1022 = -t980 * t981 + t1102;
t1023 = -t972 * t980 - t1101;
t834 = t1022 * t1108 + t1023 * t935;
t835 = -t1022 * t1105 + t1023 * t932;
t810 = t1034 * t1123 + (t1080 * t835 + t1081 * t834) * t950;
t1156 = pkin(2) * t810;
t1049 = t1156 / 0.2e1;
t941 = t980 * pkin(2);
t1052 = t810 * t941;
t804 = t807 + t810;
t1061 = t804 * t1141;
t1084 = t1007 * t950;
t1099 = t973 * t982;
t1100 = t972 * t981;
t1132 = t804 * t807;
t1147 = pkin(3) * t804;
t1183 = 0.2e1 * t982 ^ 2;
t798 = t1061 + t1156;
t864 = 0.1e1 / (pkin(2) * t981 + pkin(3) * t927 + pkin(1));
t816 = (-t973 * t993 + (-t932 * t994 + t935 * t995) * t982) * t864;
t874 = -pkin(3) + t941 + 0.2e1 * t1144;
t887 = pkin(1) - 0.2e1 * t1058;
t954 = t981 ^ 2;
t774 = (-pkin(3) * t1132 + (-0.4e1 * ((t980 * t1049 + (t953 - 0.1e1 / 0.2e1) * t1147) * t1100 + ((t1049 + t1061) * t954 - t798 / 0.2e1) * t971) * t1099 + t1190 * ((t1057 * t804 + t798 * t972) * t973 - t982 * (pkin(1) + t861) * t816) + (t1183 - 0.2e1) * t816 * (t874 * t954 + (-pkin(2) * t1102 + t887 * t980) * t981 - t1011)) * t816 / 0.2e1 + (((t1052 + (0.2e1 * t953 - 0.1e1) * t1147) * t954 - (0.2e1 * t1061 + t1156) * t971 * t1100 + t1147 - t1147 * t953) * t1183 - 0.2e1 * t816 * (t874 * t1100 + ((pkin(2) + 0.2e1 * t1141) * t954 - t909) * t971) * t1099 - 0.2e1 * t1052 - 0.2e1 * t1147 + t1190 * ((t1058 * t804 - t798 * t981) * t982 + t973 * t858 * t816)) * t810 / 0.2e1) * t1084;
t1067 = pkin(3) * t941;
t801 = t819 / 0.2e1 + t822 / 0.2e1 - t1014 / 0.2e1 + t810;
t813 = t816 ^ 2;
t783 = ((t941 + pkin(3)) * t1132 + (-t813 * ((t1174 * t953 - 0.2e1 * t1067 + t1070) * t954 - t887 * t882 + pkin(3) * t1011) + (t1004 * t804 + t1006 * t810 + 0.2e1 * t1067 * t801) * t810) * t1005) * t1084;
t786 = 0.2e1 * (-t1147 * t921 - t1156 * t972) * t864 * t816;
t849 = -t915 * t927 - t916 * t921;
t846 = -t880 * t972 - t918 * t981 + t849;
t898 = t971 * t1076;
t1033 = t864 * (t846 * t774 + t849 * t783 + (t972 * t1046 + t981 * t1074 + t894 * t1166 + t947 * t1167 + t1008 + t898 + t919 * t944 + t917 * t891 + ((t1073 + t1150) * rSges(3,2) + (-0.2e1 * t1159 + (-t926 - t980) * pkin(2)) * rSges(3,1)) * m(3)) * t786 + (0.2e1 * (-t872 * t981 - t897 * t972) * t810 + (-t1017 * t1184 - 0.2e1 * t1111 - t1120) * t816) * t810 + (-0.2e1 * t1176 * t801 - t1129) * t816 * t1172 + ((-t1117 - 0.2e1 * t1114 + (-rSges(3,1) * t921 - rSges(3,2) * t927) * t1045) * t816 + (-t915 * t921 + t916 * t927) * t804) * t804);
t1185 = t1033 * t982;
t1173 = m(3) / 0.2e1;
t1165 = t885 / 0.2e1;
t1153 = pkin(2) * t877;
t1152 = pkin(2) * t878;
t1151 = pkin(2) * t879;
t1079 = rSges(3,2) * t1184;
t1064 = rSges(3,1) * t941;
t855 = -Icges(3,3) + t898 + (-t938 - t1064) * m(3);
t890 = -m(3) * t938 - Icges(3,3);
t1087 = t774 * t855 + t783 * t890 + t786 * t849 + m(3) * t810 ^ 2 * t1153 + (t891 * t1165 + t1114 + (rSges(3,1) * t1073 + pkin(2) * t1176 + t927 * t1079 + t1153) * t1173) * t813;
t1063 = rSges(3,1) * t942;
t856 = -Icges(3,3) + t899 + (-t938 - t1063) * m(3);
t1086 = t775 * t856 + t784 * t890 + t787 * t850 + m(3) * t811 ^ 2 * t1152 + (t892 * t1165 + t1113 + (rSges(3,1) * t1072 + pkin(2) * t1177 + t929 * t1079 + t1152) * t1173) * t814;
t1062 = rSges(3,1) * t943;
t857 = -Icges(3,3) + t900 + (-t938 - t1062) * m(3);
t1085 = t776 * t857 + t785 * t890 + t788 * t851 + m(3) * t812 ^ 2 * t1151 + (t893 * t1165 + t1112 + (rSges(3,1) * t1071 + pkin(2) * t1178 + t931 * t1079 + t1151) * t1173) * t815;
t1044 = 2 * t1172;
t1024 = -(m(2) * t940) - Icges(2,3) - Icges(3,3);
t768 = t846 * t786 + (0.2e1 * t898 + (-t907 - 0.2e1 * t1064) * m(3) + t1024) * t774 + t855 * t783;
t780 = -t801 * t1044 * t1129 + (t1117 / 0.2e1 + t1114 + t1120 / 0.2e1 + t1111 + t1017 * pkin(1) + ((pkin(2) * t926 + t1159) * rSges(3,2) + (t1150 + t1162) * rSges(3,1)) * m(3)) * t813;
t1030 = (t768 + t780) * t950;
t769 = t847 * t787 + (0.2e1 * t899 + (-t907 - 0.2e1 * t1063) * m(3) + t1024) * t775 + t856 * t784;
t781 = -t802 * t1044 * t1128 + (t1116 / 0.2e1 + t1113 + t1119 / 0.2e1 + t1110 + t1016 * pkin(1) + ((pkin(2) * t928 + t1158) * rSges(3,2) + (t1149 + t1161) * rSges(3,1)) * m(3)) * t814;
t1029 = (t769 + t781) * t951;
t770 = t848 * t788 + (0.2e1 * t900 + (-t907 - 0.2e1 * t1062) * m(3) + t1024) * t776 + t857 * t785;
t782 = -t803 * t1044 * t1127 + (t1115 / 0.2e1 + t1112 + t1118 / 0.2e1 + t1109 + t1015 * pkin(1) + ((pkin(2) * t930 + t1157) * rSges(3,2) + (t1148 + t1160) * rSges(3,1)) * m(3)) * t815;
t1028 = (t770 + t782) * t952;
t1027 = t1087 * t950;
t1026 = t1086 * t951;
t1025 = t1085 * t952;
t1 = [t937 * t1187 + t936 * t1186 + t935 * t1185 + (t839 * t1028 + t837 * t1029 + t835 * t1030 + (t1025 * t842 + t1026 * t841 + t1027 * t840) * t1005) * t1007; -t934 * t1187 - t933 * t1186 - t932 * t1185 + (t838 * t1028 + t836 * t1029 + t834 * t1030 + (t1025 * t845 + t1026 * t844 + t1027 * t843) * t1005) * t1007; -t979 * t1031 - t976 * t1032 - t973 * t1033 + ((t770 / 0.2e1 + t782 / 0.2e1) * t1121 + (t769 / 0.2e1 + t781 / 0.2e1) * t1122 + (t768 / 0.2e1 + t780 / 0.2e1) * t1123 + (-t1041 * t1085 - t1042 * t1086 - t1043 * t1087) * t1005) * t1007;];
taucX  = t1;
