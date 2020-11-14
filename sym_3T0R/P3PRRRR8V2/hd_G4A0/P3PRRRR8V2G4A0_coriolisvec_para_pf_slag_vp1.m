% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V2G4A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-06 18:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G4A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:14:24
% EndTime: 2020-08-06 18:14:33
% DurationCPUTime: 9.25s
% Computational Cost: add. (53082->383), mult. (108603->740), div. (4716->10), fcn. (126546->34), ass. (0->314)
t1062 = cos(qJ(3,1));
t1040 = sin(pkin(4));
t1057 = sin(qJ(2,1));
t1111 = t1057 * t1040;
t1042 = cos(pkin(4));
t1056 = sin(qJ(3,1));
t1120 = t1042 * t1056;
t1038 = t1062 ^ 2;
t1180 = pkin(3) * t1038;
t1063 = cos(qJ(2,1));
t1070 = pkin(7) + pkin(6);
t1007 = t1063 * t1070;
t979 = pkin(2) * t1057 - t1007;
t948 = pkin(3) * t1120 + t979 * t1040;
t921 = 0.1e1 / (pkin(2) * t1120 + t948 * t1062 + t1111 * t1180);
t1060 = cos(qJ(3,2));
t1055 = sin(qJ(2,2));
t1114 = t1055 * t1040;
t1054 = sin(qJ(3,2));
t1122 = t1042 * t1054;
t1037 = t1060 ^ 2;
t1181 = pkin(3) * t1037;
t1061 = cos(qJ(2,2));
t1006 = t1061 * t1070;
t978 = pkin(2) * t1055 - t1006;
t947 = pkin(3) * t1122 + t978 * t1040;
t920 = 0.1e1 / (pkin(2) * t1122 + t947 * t1060 + t1114 * t1181);
t1058 = cos(qJ(3,3));
t1053 = sin(qJ(2,3));
t1117 = t1053 * t1040;
t1052 = sin(qJ(3,3));
t1124 = t1042 * t1052;
t1036 = t1058 ^ 2;
t1182 = pkin(3) * t1036;
t1059 = cos(qJ(2,3));
t1005 = t1059 * t1070;
t977 = pkin(2) * t1053 - t1005;
t946 = pkin(3) * t1124 + t977 * t1040;
t919 = 0.1e1 / (pkin(2) * t1124 + t946 * t1058 + t1117 * t1182);
t1083 = rSges(3,1) * t1062 - rSges(3,2) * t1056;
t1208 = t1083 * m(3);
t1084 = rSges(3,1) * t1060 - rSges(3,2) * t1054;
t1207 = t1084 * m(3);
t1085 = rSges(3,1) * t1058 - rSges(3,2) * t1052;
t1206 = t1085 * m(3);
t1043 = legFrame(3,3);
t1012 = sin(t1043);
t1018 = cos(t1043);
t1039 = sin(pkin(8));
t1041 = cos(pkin(8));
t958 = t1041 * t1012 + t1039 * t1018;
t1205 = t1040 * t958;
t1044 = legFrame(2,3);
t1013 = sin(t1044);
t1019 = cos(t1044);
t959 = t1041 * t1013 + t1039 * t1019;
t1204 = t1040 * t959;
t1045 = legFrame(1,3);
t1014 = sin(t1045);
t1020 = cos(t1045);
t960 = t1041 * t1014 + t1039 * t1020;
t1203 = t1040 * t960;
t1194 = m(3) * rSges(3,1);
t1107 = rSges(3,2) * t1194;
t1008 = -Icges(3,4) + t1107;
t1071 = pkin(2) * m(3);
t1106 = t1071 / 0.2e1;
t1099 = rSges(3,1) * t1106;
t1199 = t1008 * t1036 + t1052 * t1099;
t1198 = t1008 * t1037 + t1054 * t1099;
t1197 = t1008 * t1038 + t1056 * t1099;
t1196 = 0.2e1 * pkin(2);
t1195 = -0.2e1 * t1008;
t1072 = rSges(3,2) ^ 2;
t1073 = rSges(3,1) ^ 2;
t989 = (-t1072 + t1073) * m(3) + Icges(3,2) - Icges(3,1);
t1193 = t989 / 0.2e1;
t1064 = pkin(6) + rSges(3,3);
t1176 = t1064 * m(3);
t993 = rSges(3,2) * t1176 - Icges(3,6);
t1192 = -t993 / 0.4e1;
t994 = rSges(3,1) * t1176 - Icges(3,5);
t1191 = t994 / 0.4e1;
t986 = t1052 * rSges(3,1) + t1058 * rSges(3,2);
t1190 = m(3) * (t1085 * t1042 - t986 * t1117);
t987 = t1054 * rSges(3,1) + t1060 * rSges(3,2);
t1189 = m(3) * (t1084 * t1042 - t987 * t1114);
t988 = t1056 * rSges(3,1) + t1062 * rSges(3,2);
t1188 = m(3) * (t1083 * t1042 - t988 * t1111);
t1067 = xDP(3);
t1068 = xDP(2);
t1075 = 0.1e1 / pkin(3);
t1049 = legFrame(3,2);
t1027 = cos(t1049);
t1069 = xDP(1);
t1136 = t1027 * t1069;
t1046 = legFrame(3,1);
t1015 = sin(t1046);
t1024 = sin(t1049);
t1142 = t1015 * t1024;
t1021 = cos(t1046);
t955 = -t1039 * t1012 + t1018 * t1041;
t1151 = t955 * t1021;
t1002 = t1058 * pkin(3) + pkin(2);
t1115 = t1053 * t1070;
t1160 = t1042 * (t1002 * t1059 + t1115);
t913 = t958 * t1021 + t955 * t1142;
t970 = t1053 * t1002 - t1005;
t889 = -t913 * t1160 + (t958 * t1142 - t1151) * t970;
t1139 = t1021 * t1024;
t1166 = t1015 * t955;
t910 = -t1015 * t958 + t955 * t1139;
t892 = t910 * t1160 - (t958 * t1139 + t1166) * t970;
t904 = -t955 * t1160 + t958 * t970;
t1130 = t1040 * t1058;
t973 = t1002 * t1124;
t937 = 0.1e1 / (t970 * t1130 + t973);
t871 = (t1067 * t892 + t1068 * t889 + t904 * t1136) * t937 * t1075;
t1187 = pkin(3) * t871;
t1050 = legFrame(2,2);
t1028 = cos(t1050);
t1135 = t1028 * t1069;
t1047 = legFrame(2,1);
t1016 = sin(t1047);
t1025 = sin(t1050);
t1141 = t1016 * t1025;
t1022 = cos(t1047);
t956 = -t1039 * t1013 + t1019 * t1041;
t1150 = t956 * t1022;
t1003 = t1060 * pkin(3) + pkin(2);
t1112 = t1055 * t1070;
t1159 = t1042 * (t1003 * t1061 + t1112);
t914 = t959 * t1022 + t956 * t1141;
t971 = t1055 * t1003 - t1006;
t890 = -t914 * t1159 + (t959 * t1141 - t1150) * t971;
t1138 = t1022 * t1025;
t1165 = t1016 * t956;
t911 = -t1016 * t959 + t956 * t1138;
t893 = t911 * t1159 - (t959 * t1138 + t1165) * t971;
t905 = -t956 * t1159 + t959 * t971;
t1128 = t1040 * t1060;
t974 = t1003 * t1122;
t938 = 0.1e1 / (t971 * t1128 + t974);
t872 = (t1067 * t893 + t1068 * t890 + t905 * t1135) * t938 * t1075;
t1186 = pkin(3) * t872;
t1051 = legFrame(1,2);
t1029 = cos(t1051);
t1134 = t1029 * t1069;
t1048 = legFrame(1,1);
t1017 = sin(t1048);
t1026 = sin(t1051);
t1140 = t1017 * t1026;
t1023 = cos(t1048);
t957 = -t1039 * t1014 + t1020 * t1041;
t1149 = t957 * t1023;
t1004 = t1062 * pkin(3) + pkin(2);
t1109 = t1057 * t1070;
t1158 = t1042 * (t1004 * t1063 + t1109);
t915 = t960 * t1023 + t957 * t1140;
t972 = t1057 * t1004 - t1007;
t891 = -t915 * t1158 + (t960 * t1140 - t1149) * t972;
t1137 = t1023 * t1026;
t1164 = t1017 * t957;
t912 = -t1017 * t960 + t957 * t1137;
t894 = t912 * t1158 - (t960 * t1137 + t1164) * t972;
t906 = -t957 * t1158 + t960 * t972;
t1126 = t1040 * t1062;
t975 = t1004 * t1120;
t939 = 0.1e1 / (t972 * t1126 + t975);
t873 = (t1067 * t894 + t1068 * t891 + t906 * t1134) * t939 * t1075;
t1185 = pkin(3) * t873;
t1184 = -t1008 / 0.2e1;
t1183 = m(3) * t1042;
t1179 = t1052 * pkin(2);
t1178 = t1054 * pkin(2);
t1177 = t1056 * pkin(2);
t1076 = pkin(2) ^ 2;
t1009 = t1070 ^ 2 + t1076;
t1074 = pkin(3) ^ 2;
t1105 = t1052 * t1187;
t1172 = pkin(3) * t1196;
t1123 = t1042 * t1053;
t961 = t1039 * t1123 - t1041 * t1059;
t964 = t1039 * t1059 + t1041 * t1123;
t1082 = t1012 * t964 + t961 * t1018;
t922 = -t1012 * t961 + t964 * t1018;
t883 = (-t1015 * t1082 + t922 * t1139) * t1052 + t910 * t1130;
t886 = (-t1082 * t1021 - t922 * t1142) * t1052 - t913 * t1130;
t901 = t955 * t1130 + t1052 * (t1059 * t958 + t955 * t1123);
t1116 = t1053 * t1058;
t952 = pkin(3) * t1116 + t977;
t934 = 0.1e1 / (t952 * t1130 + t973);
t865 = -t901 * t919 * t1136 + (t1067 * t883 + t1068 * t886) * t934;
t1175 = (-t1070 * t1105 + (t1036 * t1074 + t1058 * t1172 + t1009) * t865) * t865;
t1104 = t1054 * t1186;
t1121 = t1042 * t1055;
t962 = t1039 * t1121 - t1041 * t1061;
t965 = t1039 * t1061 + t1041 * t1121;
t1081 = t1013 * t965 + t962 * t1019;
t923 = -t1013 * t962 + t965 * t1019;
t884 = (-t1016 * t1081 + t923 * t1138) * t1054 + t911 * t1128;
t887 = (-t1081 * t1022 - t923 * t1141) * t1054 - t914 * t1128;
t902 = t956 * t1128 + t1054 * (t1061 * t959 + t956 * t1121);
t1113 = t1055 * t1060;
t953 = pkin(3) * t1113 + t978;
t935 = 0.1e1 / (t953 * t1128 + t974);
t866 = -t902 * t920 * t1135 + (t1067 * t884 + t1068 * t887) * t935;
t1174 = (-t1070 * t1104 + (t1037 * t1074 + t1060 * t1172 + t1009) * t866) * t866;
t1103 = t1056 * t1185;
t1119 = t1042 * t1057;
t963 = t1039 * t1119 - t1041 * t1063;
t966 = t1039 * t1063 + t1041 * t1119;
t1080 = t1014 * t966 + t963 * t1020;
t924 = -t1014 * t963 + t966 * t1020;
t885 = (-t1017 * t1080 + t924 * t1137) * t1056 + t912 * t1126;
t888 = (-t1080 * t1023 - t924 * t1140) * t1056 - t915 * t1126;
t903 = t957 * t1126 + t1056 * (t1063 * t960 + t957 * t1119);
t1110 = t1057 * t1062;
t954 = pkin(3) * t1110 + t979;
t936 = 0.1e1 / (t954 * t1126 + t975);
t867 = -t903 * t921 * t1134 + (t1067 * t885 + t1068 * t888) * t936;
t1173 = (-t1070 * t1103 + (t1038 * t1074 + t1062 * t1172 + t1009) * t867) * t867;
t1001 = m(2) * rSges(2,1) + t1071;
t1035 = -m(1) - m(2) - m(3);
t1146 = 0.2e1 * m(3);
t1157 = t1059 * t865;
t991 = m(2) * rSges(2,2) - t1176;
t1163 = t1040 * ((t1001 + t1206) * t1059 - t1053 * t991);
t1129 = t1040 * t1059;
t1133 = t1040 * t1052;
t1154 = t865 * t1070;
t1102 = t1052 * t1154;
t856 = t1102 - t1187;
t838 = (((t1042 * t871 + t865 * t1129) * t1182 + ((-t1105 + t1154) * t1053 + pkin(2) * t1157) * t1130 + t1042 * t856) * t865 + (t871 * t1129 + (t1036 * t1042 - t1116 * t1133 - t1042) * t865) * t1187) * t919;
t1118 = t1042 * t1075;
t841 = t919 * t1118 * t1175 + (-t1042 * t1102 + (-t952 * t1133 + (pkin(2) * t1058 + t1182) * t1042) * t871) * t934 * t871;
t847 = (-t1058 * t1175 - (pkin(2) * t871 - t856 * t1058) * t1187) * t919;
t862 = t865 ^ 2;
t868 = t871 ^ 2;
t1171 = t1035 * t847 - t838 * t1163 - t841 * t1190 + ((-t862 * t1001 - (t862 + t868) * t1206) * t1053 - (t986 * t871 * t1146 + t991 * t865) * t1157) * t1040 - t868 * t986 * t1183;
t1156 = t1061 * t866;
t1162 = t1040 * ((t1001 + t1207) * t1061 - t1055 * t991);
t1127 = t1040 * t1061;
t1132 = t1040 * t1054;
t1153 = t866 * t1070;
t1101 = t1054 * t1153;
t857 = t1101 - t1186;
t839 = (((t1042 * t872 + t866 * t1127) * t1181 + ((-t1104 + t1153) * t1055 + pkin(2) * t1156) * t1128 + t1042 * t857) * t866 + (t872 * t1127 + (t1037 * t1042 - t1113 * t1132 - t1042) * t866) * t1186) * t920;
t842 = t920 * t1118 * t1174 + (-t1042 * t1101 + (-t953 * t1132 + (pkin(2) * t1060 + t1181) * t1042) * t872) * t935 * t872;
t848 = (-t1060 * t1174 - (pkin(2) * t872 - t857 * t1060) * t1186) * t920;
t863 = t866 ^ 2;
t869 = t872 ^ 2;
t1170 = t1035 * t848 - t839 * t1162 - t842 * t1189 + ((-t863 * t1001 - (t863 + t869) * t1207) * t1055 - (t987 * t872 * t1146 + t991 * t866) * t1156) * t1040 - t869 * t987 * t1183;
t1155 = t1063 * t867;
t1161 = t1040 * ((t1001 + t1208) * t1063 - t1057 * t991);
t1125 = t1040 * t1063;
t1131 = t1040 * t1056;
t1152 = t867 * t1070;
t1100 = t1056 * t1152;
t858 = t1100 - t1185;
t840 = (((t1042 * t873 + t867 * t1125) * t1180 + ((-t1103 + t1152) * t1057 + pkin(2) * t1155) * t1126 + t1042 * t858) * t867 + (t873 * t1125 + (t1038 * t1042 - t1110 * t1131 - t1042) * t867) * t1185) * t921;
t843 = t921 * t1118 * t1173 + (-t1042 * t1100 + (-t954 * t1131 + (pkin(2) * t1062 + t1180) * t1042) * t873) * t936 * t873;
t849 = (-t1062 * t1173 - (pkin(2) * t873 - t858 * t1062) * t1185) * t921;
t864 = t867 ^ 2;
t870 = t873 ^ 2;
t1169 = t1035 * t849 - t840 * t1161 - t843 * t1188 + ((-t864 * t1001 - (t864 + t870) * t1208) * t1057 - (t988 * t873 * t1146 + t991 * t867) * t1155) * t1040 - t870 * t988 * t1183;
t1148 = rSges(3,2) * t1196;
t1108 = -t1107 / 0.2e1 + Icges(3,4) / 0.2e1;
t1011 = rSges(3,2) * t1106;
t1098 = t1171 * t919;
t1097 = t1170 * t920;
t1096 = t1169 * t921;
t1031 = t1194 * t1196;
t1092 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - Icges(3,1) - Icges(2,3);
t949 = t1052 * t994 + t993 * t1058;
t992 = t1064 ^ 2 + t1072 + t1076;
t1095 = 0.4e1 * ((t1192 * t1052 + t1191 * t1058) * t871 + ((t1052 * t1193 + t1011) * t1058 + t1184 + t1199) * t865) * t871 + t847 * t1163 - (-t989 * t1036 - (t1052 * t1195 + t1031) * t1058 + (t1052 * t1148 - t992) * m(3) + t1092) * t838 - t949 * t841;
t950 = t1054 * t994 + t993 * t1060;
t1094 = 0.4e1 * ((t1192 * t1054 + t1191 * t1060) * t872 + ((t1054 * t1193 + t1011) * t1060 + t1184 + t1198) * t866) * t872 + t848 * t1162 - (-t989 * t1037 - (t1054 * t1195 + t1031) * t1060 + (t1054 * t1148 - t992) * m(3) + t1092) * t839 - t950 * t842;
t951 = t1056 * t994 + t993 * t1062;
t1093 = 0.4e1 * ((t1192 * t1056 + t1191 * t1062) * t873 + ((t1056 * t1193 + t1011) * t1062 + t1184 + t1197) * t867) * t873 + t849 * t1161 - (-t989 * t1038 - (t1056 * t1195 + t1031) * t1062 + (t1056 * t1148 - t992) * m(3) + t1092) * t840 - t951 * t843;
t1091 = t1095 * t934;
t1090 = t1094 * t935;
t1089 = t1093 * t936;
t976 = (t1073 / 0.2e1 - t1072 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t990 = -(t1072 + t1073) * m(3) - Icges(3,3);
t1088 = (-t847 * t1190 + t949 * t838 + t990 * t841 + 0.2e1 * t862 * ((t976 * t1052 + t1011) * t1058 + t1108 + t1199)) * t937;
t1087 = (-t848 * t1189 + t950 * t839 + t990 * t842 + 0.2e1 * t863 * ((t976 * t1054 + t1011) * t1060 + t1108 + t1198)) * t938;
t1086 = (-t849 * t1188 + t951 * t840 + t990 * t843 + 0.2e1 * t864 * ((t976 * t1056 + t1011) * t1062 + t1108 + t1197)) * t939;
t1079 = pkin(3) * t1133 - t1042 * t977;
t1078 = pkin(3) * t1132 - t1042 * t978;
t1077 = pkin(3) * t1131 - t1042 * t979;
t982 = pkin(2) * t1063 + t1109;
t981 = pkin(2) * t1061 + t1112;
t980 = pkin(2) * t1059 + t1115;
t933 = -t1026 * t1203 + t1042 * t1029;
t932 = -t1025 * t1204 + t1042 * t1028;
t931 = -t1024 * t1205 + t1042 * t1027;
t930 = t982 * t1039 - t1077 * t1041;
t929 = t981 * t1039 - t1078 * t1041;
t928 = t980 * t1039 - t1079 * t1041;
t927 = -t1077 * t1039 - t982 * t1041;
t926 = -t1078 * t1039 - t981 * t1041;
t925 = -t1079 * t1039 - t980 * t1041;
t909 = t1080 * t1026 + t1029 * t1111;
t908 = t1081 * t1025 + t1028 * t1114;
t907 = t1082 * t1024 + t1027 * t1117;
t900 = -t1014 * t927 + t930 * t1020;
t899 = -t1013 * t926 + t929 * t1019;
t898 = -t1012 * t925 + t928 * t1018;
t897 = t948 * t1029 + (t1014 * t930 + t1020 * t927) * t1026;
t896 = t947 * t1028 + (t1013 * t929 + t1019 * t926) * t1025;
t895 = t946 * t1027 + (t1012 * t928 + t1018 * t925) * t1024;
t1 = [(t1169 * (-((-t1063 * t957 + t960 * t1119) * t1029 - t1026 * t1111) * t1180 + ((t1077 * t960 + t982 * t957) * t1029 + t948 * t1026) * t1062 + (t1042 * t1026 + t1029 * t1203) * t1177) + t1093 * t903 * t1029) * t921 + (t1170 * (-((-t1061 * t956 + t959 * t1121) * t1028 - t1025 * t1114) * t1181 + ((t1078 * t959 + t981 * t956) * t1028 + t947 * t1025) * t1060 + (t1042 * t1025 + t1028 * t1204) * t1178) + t1094 * t902 * t1028) * t920 + (t1171 * (-((-t1059 * t955 + t958 * t1123) * t1027 - t1024 * t1117) * t1182 + ((t1079 * t958 + t980 * t955) * t1027 + t946 * t1024) * t1058 + (t1042 * t1024 + t1027 * t1205) * t1179) + t1095 * t901 * t1027) * t919 + (t904 * t1027 * t1088 + t905 * t1028 * t1087 + t906 * t1029 * t1086) * t1075; -t888 * t1089 - t887 * t1090 - t886 * t1091 + (-(t909 * t1017 - t924 * t1023) * t1180 + (-t897 * t1017 + t900 * t1023) * t1062 - (t933 * t1017 + t1040 * t1149) * t1177) * t1096 + (-(t908 * t1016 - t923 * t1022) * t1181 + (-t896 * t1016 + t899 * t1022) * t1060 - (t932 * t1016 + t1040 * t1150) * t1178) * t1097 + (-(t907 * t1015 - t922 * t1021) * t1182 + (-t895 * t1015 + t898 * t1021) * t1058 - (t931 * t1015 + t1040 * t1151) * t1179) * t1098 + (t891 * t1086 + t890 * t1087 + t889 * t1088) * t1075; -t885 * t1089 - t884 * t1090 - t883 * t1091 + ((t924 * t1017 + t909 * t1023) * t1180 + (t1017 * t900 + t897 * t1023) * t1062 + (t933 * t1023 - t1040 * t1164) * t1177) * t1096 + ((t923 * t1016 + t908 * t1022) * t1181 + (t1016 * t899 + t896 * t1022) * t1060 + (t932 * t1022 - t1040 * t1165) * t1178) * t1097 + ((t922 * t1015 + t907 * t1021) * t1182 + (t1015 * t898 + t895 * t1021) * t1058 + (t931 * t1021 - t1040 * t1166) * t1179) * t1098 + (t894 * t1086 + t893 * t1087 + t892 * t1088) * t1075;];
taucX  = t1;
