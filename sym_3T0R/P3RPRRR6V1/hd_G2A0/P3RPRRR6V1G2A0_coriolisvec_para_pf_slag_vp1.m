% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR6V1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:37
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:35:51
% EndTime: 2020-08-06 18:36:00
% DurationCPUTime: 8.68s
% Computational Cost: add. (47070->494), mult. (52671->743), div. (8139->16), fcn. (35550->106), ass. (0->371)
t1018 = cos(pkin(7));
t1030 = sin(qJ(1,1));
t1036 = cos(qJ(1,1));
t1045 = -pkin(6) - pkin(5);
t1174 = t1036 * t1045;
t1017 = sin(pkin(7));
t1175 = t1036 * t1017;
t1171 = t1045 * t1017;
t916 = -pkin(1) + t1171;
t895 = t916 * t1030;
t1280 = pkin(2) * t1175 + (t1030 * pkin(2) + t1174) * t1018 - t895;
t1028 = sin(qJ(1,2));
t1034 = cos(qJ(1,2));
t1176 = t1034 * t1045;
t1177 = t1034 * t1017;
t894 = t916 * t1028;
t1279 = pkin(2) * t1177 + (t1028 * pkin(2) + t1176) * t1018 - t894;
t1026 = sin(qJ(1,3));
t1032 = cos(qJ(1,3));
t1178 = t1032 * t1045;
t1179 = t1032 * t1017;
t893 = t916 * t1026;
t1278 = pkin(2) * t1179 + (t1026 * pkin(2) + t1178) * t1018 - t893;
t1070 = 0.1e1 / pkin(3);
t1042 = xDP(3);
t1025 = sin(qJ(3,3));
t1239 = 0.2e1 * t1045;
t1247 = -0.2e1 * pkin(2);
t1248 = -0.2e1 * pkin(1);
t1268 = 0.2e1 * pkin(2);
t976 = pkin(7) + qJ(3,3);
t944 = sin(t976);
t980 = -pkin(7) + qJ(3,3);
t946 = sin(t980);
t987 = qJ(1,3) + pkin(7);
t948 = sin(t987);
t1058 = 0.2e1 * qJ(3,3);
t994 = sin(t1058);
t1215 = (t948 * t1239 + cos(t987) * t1247 + t1032 * t1248 + (-cos(qJ(1,3) - t980) - cos(qJ(1,3) + t976)) * pkin(3)) / (t1025 * t1268 + pkin(3) * t994 + (t944 + t946) * pkin(1));
t1128 = t1042 * t1215;
t1031 = cos(qJ(3,3));
t967 = t1031 * pkin(3);
t937 = t967 + pkin(2);
t1192 = t937 * t1026;
t1010 = 0.1e1 / t1025;
t960 = t1018 * pkin(1);
t896 = 0.1e1 / (t960 + t937);
t1210 = t1010 * t896;
t1133 = ((t1178 + t1192) * t1018 - t893 + t937 * t1179) * t1210;
t1044 = xDP(1);
t1022 = legFrame(3,2);
t964 = cos(t1022);
t1187 = t964 * t1044;
t1043 = xDP(2);
t961 = sin(t1022);
t1189 = t961 * t1043;
t1277 = (t1128 / 0.6e1 - (t1187 / 0.6e1 - t1189 / 0.6e1) * t1133) * t1070;
t1069 = pkin(3) ^ 2;
t1166 = t1069 * t1070;
t1276 = (t1128 / 0.4e1 - (t1187 / 0.4e1 - t1189 / 0.4e1) * t1133) * t1166;
t1275 = (t1128 / 0.3e1 - (t1187 / 0.3e1 - t1189 / 0.3e1) * t1133) * t1070;
t1274 = (t1128 / 0.2e1 - (t1187 / 0.2e1 - t1189 / 0.2e1) * t1133) * t1070;
t1027 = sin(qJ(3,2));
t978 = pkin(7) + qJ(3,2);
t945 = sin(t978);
t982 = -pkin(7) + qJ(3,2);
t947 = sin(t982);
t990 = qJ(1,2) + pkin(7);
t949 = sin(t990);
t1061 = 0.2e1 * qJ(3,2);
t997 = sin(t1061);
t1214 = (t949 * t1239 + cos(t990) * t1247 + t1034 * t1248 + (-cos(qJ(1,2) - t982) - cos(qJ(1,2) + t978)) * pkin(3)) / (t1027 * t1268 + pkin(3) * t997 + (t945 + t947) * pkin(1));
t1127 = t1042 * t1214;
t1033 = cos(qJ(3,2));
t968 = t1033 * pkin(3);
t938 = t968 + pkin(2);
t1191 = t938 * t1028;
t1011 = 0.1e1 / t1027;
t897 = 0.1e1 / (t960 + t938);
t1209 = t1011 * t897;
t1131 = ((t1176 + t1191) * t1018 - t894 + t938 * t1177) * t1209;
t1023 = legFrame(2,2);
t965 = cos(t1023);
t1186 = t965 * t1044;
t962 = sin(t1023);
t1188 = t962 * t1043;
t1273 = (t1127 / 0.6e1 - (t1186 / 0.6e1 - t1188 / 0.6e1) * t1131) * t1070;
t1272 = (t1127 / 0.4e1 - (t1186 / 0.4e1 - t1188 / 0.4e1) * t1131) * t1166;
t1271 = (t1127 / 0.3e1 - (t1186 / 0.3e1 - t1188 / 0.3e1) * t1131) * t1070;
t1270 = (t1127 / 0.2e1 - (t1186 / 0.2e1 - t1188 / 0.2e1) * t1131) * t1070;
t1016 = t1045 ^ 2;
t1071 = pkin(2) ^ 2;
t1048 = 0.3e1 * t1071;
t1240 = t1016 + t1048;
t1072 = pkin(1) ^ 2;
t1269 = t1072 + t1240;
t1173 = t1043 / 0.2e1;
t1172 = t1044 / 0.2e1;
t1064 = 2 * qJ(3,1);
t1000 = sin(t1064);
t1009 = cos(t1064);
t1035 = cos(qJ(3,1));
t1015 = t1035 ^ 2;
t1029 = sin(qJ(3,1));
t1037 = (rSges(3,3) + pkin(5));
t1040 = -Icges(3,1) / 0.2e1;
t1049 = 0.2e1 * t1071;
t1052 = 0.2e1 * t1072;
t1065 = rSges(3,2) ^ 2;
t1066 = rSges(3,1) ^ 2;
t1161 = t1065 + t1066;
t1228 = pkin(1) * t1017;
t1074 = -Icges(3,2) / 0.2e1 + t1040 - ((2 * t1037 ^ 2) + t1049 + t1052 + t1161) * m(3) / 0.2e1 + 0.2e1 * (rSges(2,2) * m(2) - t1037 * m(3)) * t1228 - (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (-0.2e1 * rSges(2,1) * t960 - t1072) * m(2) - Icges(1,3) - Icges(2,3);
t1099 = t960 / 0.2e1 + pkin(2) / 0.2e1;
t1024 = legFrame(1,2);
t991 = qJ(1,1) + pkin(7);
t926 = t1024 + t991;
t927 = -t1024 + t991;
t880 = -sin(t926) + sin(t927);
t883 = cos(t927) + cos(t926);
t969 = t1035 * pkin(3);
t940 = t969 + pkin(2);
t898 = 0.1e1 / (t960 + t940);
t950 = sin(t991);
t846 = (-t1042 * t950 + t883 * t1172 + t880 * t1173) * t898;
t1090 = t846 * t1099;
t1096 = t1228 / 0.4e1 + t1037 / 0.4e1;
t1183 = 0.2e1 * m(3);
t1190 = t940 * t1030;
t1012 = 0.1e1 / t1029;
t1208 = t1012 * t898;
t1129 = ((t1174 + t1190) * t1018 - t895 + t940 * t1175) * t1208;
t1159 = pkin(7) + qJ(3,1);
t1160 = -pkin(7) + qJ(3,1);
t1227 = pkin(3) * t1000;
t1213 = (t950 * t1239 + cos(t991) * t1247 + t1036 * t1248 + (-cos(qJ(1,1) - t1160) - cos(qJ(1,1) + t1159)) * pkin(3)) / (t1029 * t1268 + t1227 + (sin(t1159) + sin(t1160)) * pkin(1));
t963 = sin(t1024);
t966 = cos(t1024);
t830 = (t1042 * t1213 + (t1043 * t963 - t1044 * t966) * t1129) * t1070;
t1201 = t1029 * t830;
t913 = (-t1065 + t1066) * m(3) + Icges(3,2) - Icges(3,1);
t1233 = t913 / 0.2e1;
t1234 = -t913 / 0.2e1;
t1235 = Icges(3,6) / 0.4e1;
t1236 = -Icges(3,5) / 0.4e1;
t1063 = 3 * qJ(3,1);
t1008 = cos(t1063);
t1068 = pkin(3) * t1069;
t928 = t960 + pkin(2);
t1085 = pkin(3) * (-t1035 * pkin(2) + t928 * t1008);
t1054 = 0.2e1 * pkin(7);
t971 = sin(t1054);
t1195 = t1072 * t971;
t1148 = 0.2e1 * t1195;
t1225 = t1045 * pkin(1);
t921 = t1018 * t1225;
t1089 = -0.4e1 * t921 + t1148;
t1053 = -0.2e1 * t1072;
t942 = t1072 * cos(t1054);
t1097 = -0.4e1 * t1071 + t1053 - 0.2e1 * t942 - t1069;
t1103 = -0.2e1 * t921 + t1195;
t1258 = -t1069 / 0.2e1 - t942;
t1118 = -0.2e1 * t1071 + t1258;
t1181 = 0.2e1 * pkin(3);
t1261 = -t1045 + t1228;
t1121 = t1261 * t1181;
t1137 = pkin(1) * t1171;
t939 = t969 + t1268;
t1147 = t939 * t960;
t1149 = pkin(2) * t960;
t1162 = t1016 + t1072;
t1062 = 4 * qJ(3,1);
t1167 = t1069 * cos(t1062);
t1170 = t1045 * t1072;
t1180 = 0.4e1 * pkin(3);
t970 = t1049 + t1072;
t884 = 0.4e1 * t1149 + t942 + t970;
t1193 = t884 * t1009;
t914 = pkin(2) * t1261;
t873 = t1103 + 0.2e1 * t914;
t1194 = t873 * t1009;
t1212 = t928 * sin(t1063);
t1229 = t1070 / 0.2e1;
t1238 = -0.6e1 * t1069;
t1242 = -0.2e1 * t1069 - 0.2e1 * t884;
t1246 = 0.4e1 * pkin(2);
t1141 = pkin(3) * t1201;
t820 = t1261 * t846 - t1141;
t1184 = t1072 + t1016 / 0.2e1;
t899 = 0.3e1 / 0.8e1 * t1069 + t1071 / 0.2e1 + t1184;
t918 = -0.2e1 * t1137;
t998 = sin(t1062);
t808 = ((-0.2e1 * t1069 * t1261 * t1008 + (t1261 * t939 + t1103 - t1194) * t1181) * t830 + (-0.4e1 * (0.6e1 * t1149 + t918 + t1052 - t1258 + t1240) * t1227 + t1212 * t1238 - t1068 * t998 + 0.8e1 * (-pkin(2) * t942 + t971 * t1170 - (0.3e1 / 0.4e1 * t1069 + t1048 + t1162) * t960 + (t899 - t1137) * t1247) * t1029) * t846) / (t1193 + t1167 / 0.2e1 - 0.2e1 * t1147 - t1072 + 0.2e1 * t1085 + t1118) * t846 * t1229 + (t820 * t1246 + (t1000 * t1242 - t1069 * t998 + (-t1029 * t960 - t1212) * t1180) * t830 + (-0.2e1 * t1194 + (-t1008 + t1035) * t1121 + t1089) * t846) / (0.4e1 * t1085 + t1097 - 0.4e1 * t1147 + t1167 + 0.2e1 * t1193) * t830;
t817 = (t820 - t1141) * t898 * t846;
t1112 = t1037 + t1228;
t1237 = m(3) * rSges(3,2);
t885 = -t1112 * t1237 + Icges(3,6);
t886 = t1112 * rSges(3,1) * m(3) - Icges(3,5);
t861 = t1029 * t886 - t885 * t1035;
t910 = t1035 * rSges(3,1) - t1029 * rSges(3,2);
t1158 = rSges(3,1) * t1237;
t941 = -Icges(3,4) + t1158;
t1105 = 0.4e1 * ((t1029 * t1235 + t1035 * t1236) * t830 + (t1029 * t1035 * t1233 + (t1015 - 0.1e1 / 0.2e1) * t941) * t846 + ((t1035 * t1090 - t1096 * t1201) * rSges(3,2) + (t1096 * t1035 * t830 + t1029 * t1090) * rSges(3,1)) * m(3)) * t830 - (t1009 * t1234 + t941 * t1000 + (-t910 * t960 + (-t960 - t910) * pkin(2)) * t1183 + t1074) * t817 - t861 * t808;
t1267 = -t1105 * t898 / 0.2e1;
t1003 = cos(t1058);
t1013 = t1031 ^ 2;
t1126 = t948 * t896 * t1042;
t922 = t1022 + t987;
t923 = -t1022 + t987;
t878 = -sin(t922) + sin(t923);
t866 = t878 * t896 * t1173;
t881 = cos(t923) + cos(t922);
t868 = t881 * t896 * t1172;
t844 = t866 + t868 - t1126;
t1092 = t844 * t1099;
t828 = (t1128 + (-t1187 + t1189) * t1133) * t1070;
t1207 = t1025 * t828;
t1245 = 0.4e1 * t828;
t1057 = 0.3e1 * qJ(3,3);
t1002 = cos(t1057);
t1076 = (0.3e1 / 0.4e1 * t868 + 0.3e1 / 0.4e1 * t866 - 0.3e1 / 0.4e1 * t1126) * t1069 + t1269 * t844;
t1098 = -0.4e1 * t1149 + t1118;
t1115 = pkin(3) * t1148;
t1226 = pkin(3) * t1045;
t1122 = t1226 * t1246;
t1150 = pkin(3) * t1225;
t1123 = 0.4e1 * t1150;
t1124 = -0.4e1 * t1150;
t1135 = -t1226 / 0.2e1;
t1136 = pkin(2) * t1181;
t1138 = -0.12e2 * pkin(1) * pkin(2) * pkin(3);
t1139 = t1069 * t1239;
t1152 = -0.4e1 * t960;
t1153 = pkin(3) * t1053;
t1154 = -0.16e2 * pkin(2) * t899;
t1155 = -0.4e1 * pkin(2) * t1072;
t1156 = pkin(2) * t1238;
t1157 = -0.3e1 * pkin(1) * t1069;
t1164 = t1002 - t1031;
t1165 = -0.8e1 * pkin(3) * (t1069 / 0.4e1 + 0.3e1 / 0.2e1 * t1071 + t1184);
t1056 = 0.4e1 * qJ(3,3);
t1169 = t1069 * cos(t1056);
t1211 = t914 - t921;
t1241 = 0.2e1 * t884;
t1243 = -0.2e1 * t873;
t1249 = -0.4e1 * pkin(1);
t1259 = 0.4e1 * t1170;
t1260 = 0.8e1 * pkin(2) * t1225;
t1143 = pkin(3) * t1207;
t818 = t1261 * t844 - t1143;
t821 = t844 - t1274;
t822 = t844 + t1274;
t932 = 0.2e1 * t976;
t934 = 0.2e1 * t980;
t975 = pkin(7) + t1058;
t951 = cos(t975);
t952 = cos(t976);
t979 = -pkin(7) + t1058;
t955 = cos(t979);
t956 = cos(t980);
t973 = t1054 + qJ(3,3);
t1055 = -0.2e1 * pkin(7);
t983 = t1055 + qJ(3,3);
t985 = t1057 + pkin(7);
t986 = t1057 - pkin(7);
t992 = sin(t1056);
t993 = sin(t1057);
t803 = (t822 * t955 * t1123 + t821 * t951 * t1124 + pkin(3) * (t1031 * t1135 + t1211) * t1245 + ((t1076 + t1276) * t946 + (t1076 - t1276) * t944) * t1249 + (t1025 * t1154 - t1068 * t992 + t993 * t1156 + t994 * t1165) * t844 + ((t844 - t1275) * sin(t986) + (t844 + t1275) * sin(t985)) * t1157 + (t821 * sin(t934) + t822 * sin(t932)) * t1153 + ((t844 - t1277) * sin(t979) + (t844 + t1277) * sin(t975)) * t1138 + (cos(t983) - cos(t973)) * t844 * t1259 + (sin(t983) + sin(t973)) * t844 * t1155 + (-t952 + t956) * t844 * t1260 + (t1002 * t1139 + t1003 * t1122 + t1115) * t828) / (t970 * t1003 + t1169 / 0.2e1 + (cos(t934) / 0.2e1 + cos(t932) / 0.2e1 - 0.1e1) * t1072 + t1164 * t1136 + ((t951 + t955) * t1268 + (cos(t986) + cos(t985) - t952 - t956) * pkin(3)) * pkin(1) + t1098) * t844 * t1229 + (t818 * t1246 + (-t1069 * t992 + t994 * t1242 + (-t1025 * t960 - t928 * t993) * t1180) * t828 + (t1003 * t1243 - t1164 * t1121 + t1089) * t844) / (t1003 * t1241 + t1169 + (t967 + t1268) * t1152 + (-t1031 * pkin(2) + t928 * t1002) * t1180 + t1097) * t828;
t815 = (t818 - t1143) * t896 * t844;
t859 = t1025 * t886 - t885 * t1031;
t906 = t1031 * rSges(3,1) - t1025 * rSges(3,2);
t1108 = ((t1025 * t1235 + t1031 * t1236) * t828 + (t1025 * t1031 * t1233 + (t1013 - 0.1e1 / 0.2e1) * t941) * t844 + ((t1031 * t1092 - t1096 * t1207) * rSges(3,2) + (t1096 * t1031 * t828 + t1025 * t1092) * rSges(3,1)) * m(3)) * t1245 - (t1003 * t1234 + t941 * t994 + (-t906 * t960 + (-t960 - t906) * pkin(2)) * t1183 + t1074) * t815 - t859 * t803;
t1266 = -t1108 * t896 / 0.2e1;
t1006 = cos(t1061);
t1014 = t1033 ^ 2;
t1125 = t949 * t897 * t1042;
t924 = t1023 + t990;
t925 = -t1023 + t990;
t879 = -sin(t924) + sin(t925);
t867 = t879 * t897 * t1173;
t882 = cos(t925) + cos(t924);
t869 = t882 * t897 * t1172;
t845 = t867 + t869 - t1125;
t1091 = t845 * t1099;
t829 = (t1127 + (-t1186 + t1188) * t1131) * t1070;
t1204 = t1027 * t829;
t1244 = 0.4e1 * t829;
t1060 = 0.3e1 * qJ(3,2);
t1005 = cos(t1060);
t1075 = (0.3e1 / 0.4e1 * t869 + 0.3e1 / 0.4e1 * t867 - 0.3e1 / 0.4e1 * t1125) * t1069 + t1269 * t845;
t1163 = t1005 - t1033;
t1059 = 0.4e1 * qJ(3,2);
t1168 = t1069 * cos(t1059);
t1142 = pkin(3) * t1204;
t819 = t1261 * t845 - t1142;
t823 = t845 - t1270;
t824 = t845 + t1270;
t933 = 0.2e1 * t978;
t935 = 0.2e1 * t982;
t977 = pkin(7) + t1061;
t953 = cos(t977);
t954 = cos(t978);
t981 = -pkin(7) + t1061;
t957 = cos(t981);
t958 = cos(t982);
t974 = t1054 + qJ(3,2);
t984 = t1055 + qJ(3,2);
t988 = t1060 + pkin(7);
t989 = t1060 - pkin(7);
t995 = sin(t1059);
t996 = sin(t1060);
t804 = (t824 * t957 * t1123 + t823 * t953 * t1124 + pkin(3) * (t1033 * t1135 + t1211) * t1244 + ((t1075 + t1272) * t947 + (t1075 - t1272) * t945) * t1249 + (t1027 * t1154 - t1068 * t995 + t996 * t1156 + t997 * t1165) * t845 + ((t845 - t1271) * sin(t989) + (t845 + t1271) * sin(t988)) * t1157 + (t823 * sin(t935) + t824 * sin(t933)) * t1153 + ((t845 - t1273) * sin(t981) + (t845 + t1273) * sin(t977)) * t1138 + (cos(t984) - cos(t974)) * t845 * t1259 + (sin(t984) + sin(t974)) * t845 * t1155 + (-t954 + t958) * t845 * t1260 + (t1005 * t1139 + t1006 * t1122 + t1115) * t829) / (t970 * t1006 + t1168 / 0.2e1 + (cos(t935) / 0.2e1 + cos(t933) / 0.2e1 - 0.1e1) * t1072 + t1163 * t1136 + ((t953 + t957) * t1268 + (cos(t989) + cos(t988) - t954 - t958) * pkin(3)) * pkin(1) + t1098) * t845 * t1229 + (t819 * t1246 + (-t1069 * t995 + t997 * t1242 + (-t1027 * t960 - t928 * t996) * t1180) * t829 + (t1006 * t1243 - t1163 * t1121 + t1089) * t845) / (t1006 * t1241 + t1168 + (t968 + t1268) * t1152 + (-t1033 * pkin(2) + t928 * t1005) * t1180 + t1097) * t829;
t816 = (t819 - t1142) * t897 * t845;
t860 = t1027 * t886 - t885 * t1033;
t908 = t1033 * rSges(3,1) - t1027 * rSges(3,2);
t1107 = ((t1027 * t1235 + t1033 * t1236) * t829 + (t1027 * t1033 * t1233 + (t1014 - 0.1e1 / 0.2e1) * t941) * t845 + ((t1033 * t1091 - t1096 * t1204) * rSges(3,2) + (t1096 * t1033 * t829 + t1027 * t1091) * rSges(3,1)) * m(3)) * t1244 - (t1006 * t1234 + t941 * t997 + (-t908 * t960 + (-t960 - t908) * pkin(2)) * t1183 + t1074) * t816 - t860 * t804;
t1265 = -t897 * t1107 / 0.2e1;
t1140 = m(3) * t928 / 0.2e1;
t1104 = rSges(3,1) * t1140;
t1185 = Icges(3,4) / 0.2e1 - t1158 / 0.2e1;
t1230 = m(3) * t910;
t1196 = t1035 * t846;
t827 = t830 ^ 2;
t876 = t1071 + t918 + 0.2e1 * t1149 + t1162;
t814 = (pkin(3) * (-t928 - t969) * t827 + (-(t1015 * t1069 + t876) * t846 + (-t928 * t1196 + t1201 * t1261) * t1181) * t1196) * t1208;
t892 = (t1066 / 0.2e1 - t1065 / 0.2e1) * m(3) + t1040 + Icges(3,2) / 0.2e1;
t904 = rSges(3,2) * t1140;
t915 = -t1161 * m(3) - Icges(3,3);
t1106 = 0.2e1 * t846 ^ 2 * (t941 * t1015 + (t892 * t1029 + t904) * t1035 + t1029 * t1104 + t1185) - t814 * t1230 + t915 * t808 + t861 * t817;
t1264 = t1106 * t1129;
t1231 = m(3) * t908;
t1197 = t1033 * t845;
t826 = t829 ^ 2;
t813 = (pkin(3) * (-t928 - t968) * t826 + (-(t1014 * t1069 + t876) * t845 + (-t928 * t1197 + t1204 * t1261) * t1181) * t1197) * t1209;
t1109 = 0.2e1 * t845 ^ 2 * (t941 * t1014 + (t892 * t1027 + t904) * t1033 + t1027 * t1104 + t1185) - t813 * t1231 + t915 * t804 + t860 * t816;
t1263 = t1109 * t1131;
t1232 = m(3) * t906;
t1198 = t1031 * t844;
t825 = t828 ^ 2;
t812 = (pkin(3) * (-t928 - t967) * t825 + (-(t1013 * t1069 + t876) * t844 + (-t928 * t1198 + t1207 * t1261) * t1181) * t1198) * t1210;
t1110 = 0.2e1 * t844 ^ 2 * (t941 * t1013 + (t892 * t1025 + t904) * t1031 + t1025 * t1104 + t1185) - t812 * t1232 + t915 * t803 + t859 * t815;
t1262 = t1110 * t1133;
t1218 = t825 * (t1025 * rSges(3,1) + t1031 * rSges(3,2));
t1217 = t826 * (t1027 * rSges(3,1) + t1033 * rSges(3,2));
t1216 = t827 * (t1029 * rSges(3,1) + t1035 * rSges(3,2));
t1206 = t1025 * t961;
t1205 = t1025 * t964;
t1203 = t1027 * t962;
t1202 = t1027 * t965;
t1200 = t1029 * t963;
t1199 = t1029 * t966;
t1146 = pkin(3) * t1013 * (t1026 * t1018 + t1179);
t1145 = pkin(3) * t1014 * (t1028 * t1018 + t1177);
t1144 = pkin(3) * t1015 * (t1030 * t1018 + t1175);
t1047 = -m(2) - m(3);
t799 = t1047 * t812 - t803 * t1232;
t1134 = t799 * t1210;
t800 = t1047 * t813 - t804 * t1231;
t1132 = t800 * t1209;
t806 = t1047 * t814 - t808 * t1230;
t1130 = t806 * t1208;
t1102 = t1210 * t1218;
t1101 = t1209 * t1217;
t1100 = t1208 * t1216;
t838 = t966 * t1144 + (pkin(3) * t1200 + t1280 * t966) * t1035 + t928 * t1200;
t837 = -t963 * t1144 + (pkin(3) * t1199 - t1280 * t963) * t1035 + t928 * t1199;
t836 = t965 * t1145 + (pkin(3) * t1203 + t1279 * t965) * t1033 + t928 * t1203;
t835 = -t962 * t1145 + (pkin(3) * t1202 - t1279 * t962) * t1033 + t928 * t1202;
t834 = t964 * t1146 + (pkin(3) * t1206 + t1278 * t964) * t1031 + t928 * t1206;
t833 = -t961 * t1146 + (pkin(3) * t1205 - t1278 * t961) * t1031 + t928 * t1205;
t1 = [t834 * t1134 + t836 * t1132 + t838 * t1130 + (-t838 * t1100 - t836 * t1101 - t834 * t1102) * m(3) + (-t964 * t1262 - t965 * t1263 - t966 * t1264) * t1070 + t883 * t1267 + t881 * t1266 + t882 * t1265; t833 * t1134 + t835 * t1132 + t837 * t1130 + (-t1100 * t837 - t1101 * t835 - t1102 * t833) * m(3) + (t961 * t1262 + t962 * t1263 + t963 * t1264) * t1070 + t880 * t1267 + t878 * t1266 + t879 * t1265; (t1105 * t950 + (-m(3) * t1216 + t806) * ((-t1030 * t1045 + t940 * t1036) * t1018 - t916 * t1036 - t1017 * t1190) * t1035 * t1012) * t898 + (t1107 * t949 + (-m(3) * t1217 + t800) * ((-t1028 * t1045 + t938 * t1034) * t1018 - t916 * t1034 - t1017 * t1191) * t1033 * t1011) * t897 + (t1108 * t948 + (-m(3) * t1218 + t799) * ((-t1026 * t1045 + t937 * t1032) * t1018 - t916 * t1032 - t1017 * t1192) * t1031 * t1010) * t896 + (t1106 * t1213 + t1109 * t1214 + t1110 * t1215) * t1070;];
taucX  = t1;
