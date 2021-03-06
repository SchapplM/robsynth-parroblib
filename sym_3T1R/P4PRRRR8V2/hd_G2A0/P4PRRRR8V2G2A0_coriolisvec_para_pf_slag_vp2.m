% Calculate vector of centrifugal and coriolis load on the joints for
% P4PRRRR8V2G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [4x3]
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
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [4x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:22
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:17:00
% EndTime: 2020-08-07 11:17:18
% DurationCPUTime: 19.15s
% Computational Cost: add. (75040->609), mult. (162680->1174), div. (9024->9), fcn. (155268->30), ass. (0->443)
t1156 = cos(pkin(4));
t1176 = cos(qJ(3,1));
t1170 = sin(qJ(3,1));
t1405 = mrSges(3,1) * t1170;
t1230 = mrSges(3,2) * t1176 + t1405;
t1246 = t1176 * mrSges(3,1) - mrSges(3,2) * t1170;
t1154 = sin(pkin(4));
t1171 = sin(qJ(2,1));
t1318 = t1154 * t1171;
t1046 = t1246 * t1156 - t1230 * t1318;
t1126 = mrSges(3,2) * pkin(6) - Ifges(3,6);
t1127 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t1083 = -t1126 * t1176 - t1170 * t1127;
t1151 = t1176 ^ 2;
t1161 = Ifges(3,1) - Ifges(3,2);
t1178 = pkin(2) * mrSges(3,2);
t1427 = -2 * Ifges(3,4);
t1431 = Ifges(3,4) + t1151 * t1427 + (-t1161 * t1170 + t1178) * t1176;
t1177 = cos(qJ(2,1));
t1184 = pkin(7) + pkin(6);
t1099 = pkin(2) * t1171 - t1177 * t1184;
t1299 = t1156 * t1170;
t1058 = pkin(3) * t1299 + t1099 * t1154;
t1413 = pkin(3) * t1151;
t1026 = 0.1e1 / (pkin(2) * t1299 + t1058 * t1176 + t1318 * t1413);
t1102 = pkin(2) * t1177 + t1171 * t1184;
t1153 = sin(pkin(8));
t1155 = cos(pkin(8));
t1295 = t1156 * t1177;
t1311 = t1155 * t1156;
t1409 = pkin(3) * t1176;
t1019 = (t1153 * t1171 - t1155 * t1295) * t1409 - t1102 * t1311 + t1099 * t1153;
t1329 = t1153 * t1156;
t1022 = (t1153 * t1295 + t1155 * t1171) * t1409 + t1102 * t1329 + t1099 * t1155;
t1181 = xDP(3);
t1186 = xP(4);
t1144 = sin(t1186);
t1145 = cos(t1186);
t1192 = koppelP(1,2);
t1196 = koppelP(1,1);
t1090 = t1144 * t1196 + t1145 * t1192;
t1094 = -t1144 * t1192 + t1145 * t1196;
t1165 = legFrame(1,2);
t1135 = sin(t1165);
t1139 = cos(t1165);
t1180 = xDP(4);
t1182 = xDP(2);
t1183 = xDP(1);
t1222 = (t1094 * t1180 + t1182) * t1135 - (-t1090 * t1180 + t1183) * t1139;
t1199 = 0.1e1 / pkin(3);
t1350 = t1026 * t1199;
t983 = (t1019 * t1181 + t1222 * t1022) * t1350;
t1423 = pkin(3) * t983;
t1274 = t1170 * t1423;
t1286 = t1171 * t1176;
t1312 = t1154 * t1177;
t1313 = t1154 * t1176;
t1319 = t1154 * t1170;
t1298 = t1156 * t1171;
t1076 = t1153 * t1298 - t1155 * t1177;
t1040 = t1076 * t1170 + t1153 * t1313;
t1079 = t1153 * t1177 + t1155 * t1298;
t1307 = t1155 * t1176;
t1043 = -t1079 * t1170 - t1154 * t1307;
t994 = (t1222 * t1040 + t1043 * t1181) * t1026;
t1375 = t1184 * t994;
t1379 = t1177 * t994;
t1387 = t1026 * t983;
t1258 = t1170 * t1375;
t935 = t1258 - t1423;
t898 = ((t1156 * t983 + t994 * t1312) * t1413 + ((-t1274 + t1375) * t1171 + pkin(2) * t1379) * t1313 + t935 * t1156) * t1026 * t994 + pkin(3) * (t983 * t1312 + (t1151 * t1156 - t1286 * t1319 - t1156) * t994) * t1387;
t1070 = pkin(3) * t1286 + t1099;
t1294 = t1156 * t1199;
t1200 = pkin(2) ^ 2;
t1123 = t1184 ^ 2 + t1200;
t1198 = pkin(3) ^ 2;
t1395 = 0.2e1 * pkin(2) * pkin(3);
t1396 = (-t1184 * t1274 + (t1151 * t1198 + t1176 * t1395 + t1123) * t994) * t994;
t902 = t1026 * t1294 * t1396 + (-t1156 * t1258 + (-t1070 * t1319 + (pkin(2) * t1176 + t1413) * t1156) * t983) / (t1070 * t1313 + (pkin(2) + t1409) * t1299) * t983;
t910 = (-t1176 * t1396 - (pkin(2) * t983 - t1176 * t935) * t1423) * t1026;
t991 = t994 ^ 2;
t1469 = -Ifges(3,3) * t902 - t1046 * t910 - t1083 * t898 + (pkin(2) * t1405 + t1431) * t991;
t1174 = cos(qJ(3,2));
t1168 = sin(qJ(3,2));
t1406 = mrSges(3,1) * t1168;
t1231 = mrSges(3,2) * t1174 + t1406;
t1247 = t1174 * mrSges(3,1) - mrSges(3,2) * t1168;
t1169 = sin(qJ(2,2));
t1320 = t1154 * t1169;
t1045 = t1247 * t1156 - t1231 * t1320;
t1082 = -t1126 * t1174 - t1168 * t1127;
t1150 = t1174 ^ 2;
t1430 = Ifges(3,4) + t1150 * t1427 + (-t1161 * t1168 + t1178) * t1174;
t1175 = cos(qJ(2,2));
t1098 = pkin(2) * t1169 - t1175 * t1184;
t1301 = t1156 * t1168;
t1057 = pkin(3) * t1301 + t1098 * t1154;
t1414 = pkin(3) * t1150;
t1025 = 0.1e1 / (pkin(2) * t1301 + t1057 * t1174 + t1320 * t1414);
t1101 = pkin(2) * t1175 + t1169 * t1184;
t1296 = t1156 * t1175;
t1410 = pkin(3) * t1174;
t1018 = (t1153 * t1169 - t1155 * t1296) * t1410 - t1101 * t1311 + t1098 * t1153;
t1021 = (t1153 * t1296 + t1155 * t1169) * t1410 + t1101 * t1329 + t1098 * t1155;
t1191 = koppelP(2,2);
t1195 = koppelP(2,1);
t1089 = t1144 * t1195 + t1145 * t1191;
t1093 = -t1144 * t1191 + t1145 * t1195;
t1164 = legFrame(2,2);
t1134 = sin(t1164);
t1138 = cos(t1164);
t1223 = (t1093 * t1180 + t1182) * t1134 - (-t1089 * t1180 + t1183) * t1138;
t1352 = t1025 * t1199;
t982 = (t1018 * t1181 + t1223 * t1021) * t1352;
t1424 = pkin(3) * t982;
t1275 = t1168 * t1424;
t1287 = t1169 * t1174;
t1314 = t1154 * t1175;
t1315 = t1154 * t1174;
t1321 = t1154 * t1168;
t1300 = t1156 * t1169;
t1075 = t1153 * t1300 - t1155 * t1175;
t1039 = t1075 * t1168 + t1153 * t1315;
t1078 = t1153 * t1175 + t1155 * t1300;
t1308 = t1155 * t1174;
t1042 = -t1078 * t1168 - t1154 * t1308;
t993 = (t1223 * t1039 + t1042 * t1181) * t1025;
t1376 = t1184 * t993;
t1380 = t1175 * t993;
t1389 = t1025 * t982;
t1259 = t1168 * t1376;
t934 = t1259 - t1424;
t897 = ((t1156 * t982 + t993 * t1314) * t1414 + ((-t1275 + t1376) * t1169 + pkin(2) * t1380) * t1315 + t934 * t1156) * t1025 * t993 + pkin(3) * (t982 * t1314 + (t1150 * t1156 - t1287 * t1321 - t1156) * t993) * t1389;
t1069 = pkin(3) * t1287 + t1098;
t1397 = (-t1184 * t1275 + (t1150 * t1198 + t1174 * t1395 + t1123) * t993) * t993;
t901 = t1025 * t1294 * t1397 + (-t1156 * t1259 + (-t1069 * t1321 + (pkin(2) * t1174 + t1414) * t1156) * t982) / (t1069 * t1315 + (pkin(2) + t1410) * t1301) * t982;
t909 = (-t1174 * t1397 - (pkin(2) * t982 - t1174 * t934) * t1424) * t1025;
t990 = t993 ^ 2;
t1468 = -Ifges(3,3) * t901 - t1045 * t909 - t1082 * t897 + (pkin(2) * t1406 + t1430) * t990;
t1159 = cos(qJ(3,4));
t1157 = sin(qJ(3,4));
t1408 = mrSges(3,1) * t1157;
t1233 = mrSges(3,2) * t1159 + t1408;
t1249 = t1159 * mrSges(3,1) - mrSges(3,2) * t1157;
t1158 = sin(qJ(2,4));
t1326 = t1154 * t1158;
t1035 = t1249 * t1156 - t1233 * t1326;
t1073 = -t1126 * t1159 - t1157 * t1127;
t1146 = t1159 ^ 2;
t1428 = Ifges(3,4) + t1146 * t1427 + (-t1157 * t1161 + t1178) * t1159;
t1160 = cos(qJ(2,4));
t1095 = pkin(2) * t1158 - t1160 * t1184;
t1306 = t1156 * t1157;
t1055 = pkin(3) * t1306 + t1095 * t1154;
t1416 = pkin(3) * t1146;
t1023 = 0.1e1 / (pkin(2) * t1306 + t1055 * t1159 + t1326 * t1416);
t1096 = pkin(2) * t1160 + t1158 * t1184;
t1304 = t1156 * t1160;
t1412 = pkin(3) * t1159;
t1015 = (t1153 * t1158 - t1155 * t1304) * t1412 - t1096 * t1311 + t1095 * t1153;
t1016 = (t1153 * t1304 + t1155 * t1158) * t1412 + t1096 * t1329 + t1095 * t1155;
t1189 = koppelP(4,2);
t1193 = koppelP(4,1);
t1087 = t1144 * t1193 + t1145 * t1189;
t1091 = -t1144 * t1189 + t1145 * t1193;
t1162 = legFrame(4,2);
t1132 = sin(t1162);
t1136 = cos(t1162);
t1225 = (t1091 * t1180 + t1182) * t1132 - (-t1087 * t1180 + t1183) * t1136;
t1356 = t1023 * t1199;
t976 = (t1015 * t1181 + t1225 * t1016) * t1356;
t1426 = pkin(3) * t976;
t1277 = t1157 * t1426;
t1292 = t1158 * t1159;
t1324 = t1154 * t1160;
t1325 = t1154 * t1159;
t1327 = t1154 * t1157;
t1305 = t1156 * t1158;
t1071 = t1153 * t1305 - t1155 * t1160;
t1036 = t1071 * t1157 + t1153 * t1325;
t1072 = t1153 * t1160 + t1155 * t1305;
t1310 = t1155 * t1159;
t1037 = -t1072 * t1157 - t1154 * t1310;
t988 = (t1225 * t1036 + t1037 * t1181) * t1023;
t1378 = t1184 * t988;
t1382 = t1160 * t988;
t1393 = t1023 * t976;
t1261 = t1157 * t1378;
t931 = t1261 - t1426;
t895 = ((t1156 * t976 + t988 * t1324) * t1416 + ((-t1277 + t1378) * t1158 + pkin(2) * t1382) * t1325 + t931 * t1156) * t1023 * t988 + pkin(3) * (t976 * t1324 + (t1146 * t1156 - t1292 * t1327 - t1156) * t988) * t1393;
t1067 = pkin(3) * t1292 + t1095;
t1399 = (-t1184 * t1277 + (t1146 * t1198 + t1159 * t1395 + t1123) * t988) * t988;
t899 = t1023 * t1294 * t1399 + (-t1156 * t1261 + (-t1067 * t1327 + (pkin(2) * t1159 + t1416) * t1156) * t976) / (t1067 * t1325 + (pkin(2) + t1412) * t1306) * t976;
t907 = (-t1159 * t1399 - (pkin(2) * t976 - t1159 * t931) * t1426) * t1023;
t987 = t988 ^ 2;
t1467 = -Ifges(3,3) * t899 - t1035 * t907 - t1073 * t895 + (pkin(2) * t1408 + t1428) * t987;
t1172 = cos(qJ(3,3));
t1166 = sin(qJ(3,3));
t1407 = mrSges(3,1) * t1166;
t1232 = mrSges(3,2) * t1172 + t1407;
t1248 = t1172 * mrSges(3,1) - mrSges(3,2) * t1166;
t1167 = sin(qJ(2,3));
t1322 = t1154 * t1167;
t1044 = t1248 * t1156 - t1232 * t1322;
t1081 = -t1126 * t1172 - t1166 * t1127;
t1149 = t1172 ^ 2;
t1429 = Ifges(3,4) + t1149 * t1427 + (-t1161 * t1166 + t1178) * t1172;
t1173 = cos(qJ(2,3));
t1097 = pkin(2) * t1167 - t1173 * t1184;
t1303 = t1156 * t1166;
t1056 = pkin(3) * t1303 + t1097 * t1154;
t1415 = pkin(3) * t1149;
t1024 = 0.1e1 / (pkin(2) * t1303 + t1056 * t1172 + t1322 * t1415);
t1100 = pkin(2) * t1173 + t1167 * t1184;
t1297 = t1156 * t1173;
t1411 = pkin(3) * t1172;
t1017 = (t1153 * t1167 - t1155 * t1297) * t1411 - t1100 * t1311 + t1097 * t1153;
t1020 = (t1153 * t1297 + t1155 * t1167) * t1411 + t1100 * t1329 + t1097 * t1155;
t1190 = koppelP(3,2);
t1194 = koppelP(3,1);
t1088 = t1144 * t1194 + t1145 * t1190;
t1092 = -t1144 * t1190 + t1145 * t1194;
t1163 = legFrame(3,2);
t1133 = sin(t1163);
t1137 = cos(t1163);
t1224 = (t1092 * t1180 + t1182) * t1133 - (-t1088 * t1180 + t1183) * t1137;
t1354 = t1024 * t1199;
t981 = (t1017 * t1181 + t1224 * t1020) * t1354;
t1425 = pkin(3) * t981;
t1276 = t1166 * t1425;
t1288 = t1167 * t1172;
t1316 = t1154 * t1173;
t1317 = t1154 * t1172;
t1323 = t1154 * t1166;
t1302 = t1156 * t1167;
t1074 = t1153 * t1302 - t1155 * t1173;
t1038 = t1074 * t1166 + t1153 * t1317;
t1077 = t1153 * t1173 + t1155 * t1302;
t1309 = t1155 * t1172;
t1041 = -t1077 * t1166 - t1154 * t1309;
t992 = (t1224 * t1038 + t1041 * t1181) * t1024;
t1377 = t1184 * t992;
t1381 = t1173 * t992;
t1391 = t1024 * t981;
t1260 = t1166 * t1377;
t933 = t1260 - t1425;
t896 = ((t1156 * t981 + t992 * t1316) * t1415 + ((-t1276 + t1377) * t1167 + pkin(2) * t1381) * t1317 + t933 * t1156) * t1024 * t992 + pkin(3) * (t981 * t1316 + (t1149 * t1156 - t1288 * t1323 - t1156) * t992) * t1391;
t1068 = pkin(3) * t1288 + t1097;
t1398 = (-t1184 * t1276 + (t1149 * t1198 + t1172 * t1395 + t1123) * t992) * t992;
t900 = t1024 * t1294 * t1398 + (-t1156 * t1260 + (-t1068 * t1323 + (pkin(2) * t1172 + t1415) * t1156) * t981) / (t1068 * t1317 + (pkin(2) + t1411) * t1303) * t981;
t908 = (-t1172 * t1398 - (pkin(2) * t981 - t1172 * t933) * t1425) * t1024;
t989 = t992 ^ 2;
t1466 = -Ifges(3,3) * t900 - t1044 * t908 - t1081 * t896 + (pkin(2) * t1407 + t1429) * t989;
t1465 = t1469 * t1350;
t1464 = t1468 * t1352;
t1463 = t1466 * t1354;
t1462 = t1467 * t1356;
t1179 = mrSges(3,1) * pkin(2);
t1205 = Ifges(3,1) + Ifges(2,3) + (pkin(6) ^ 2 + t1200) * m(3) + 0.2e1 * pkin(6) * mrSges(3,3);
t1374 = -0.2e1 * t1178;
t1032 = -t1161 * t1149 + 0.2e1 * (Ifges(3,4) * t1166 + t1179) * t1172 + t1166 * t1374 + t1205;
t1124 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t1140 = m(3) * pkin(2) + mrSges(2,1);
t1336 = ((t1248 + t1140) * t1173 + t1167 * t1124) * t1154;
t886 = -t1032 * t896 - t1081 * t900 - t908 * t1336;
t1421 = t1127 / 0.2e1;
t1422 = -t1126 / 0.2e1;
t912 = (t1422 * t1166 + t1421 * t1172) * t981 + (t1179 * t1166 + t1429) * t992;
t1461 = -t1024 * t886 + 0.2e1 * t912 * t1391;
t1031 = -t1161 * t1146 + 0.2e1 * (Ifges(3,4) * t1157 + t1179) * t1159 + t1157 * t1374 + t1205;
t1337 = ((t1249 + t1140) * t1160 + t1158 * t1124) * t1154;
t883 = -t1031 * t895 - t1073 * t899 - t907 * t1337;
t911 = (t1422 * t1157 + t1421 * t1159) * t976 + (t1179 * t1157 + t1428) * t988;
t1460 = -t1023 * t883 + 0.2e1 * t911 * t1393;
t1033 = -t1161 * t1150 + 0.2e1 * (Ifges(3,4) * t1168 + t1179) * t1174 + t1168 * t1374 + t1205;
t1335 = ((t1247 + t1140) * t1175 + t1169 * t1124) * t1154;
t887 = -t1033 * t897 - t1082 * t901 - t909 * t1335;
t913 = (t1422 * t1168 + t1421 * t1174) * t982 + (t1179 * t1168 + t1430) * t993;
t1459 = t1025 * t887 - 0.2e1 * t913 * t1389;
t1034 = -t1161 * t1151 + 0.2e1 * (Ifges(3,4) * t1170 + t1179) * t1176 + t1170 * t1374 + t1205;
t1334 = ((t1246 + t1140) * t1177 + t1171 * t1124) * t1154;
t888 = -t1034 * t898 - t1083 * t902 - t910 * t1334;
t914 = (t1422 * t1170 + t1421 * t1176) * t983 + (t1179 * t1170 + t1431) * t994;
t1458 = t1026 * t888 - 0.2e1 * t914 * t1387;
t1451 = -0.2e1 * t976;
t1450 = -0.2e1 * t981;
t1449 = -0.2e1 * t982;
t1448 = -0.2e1 * t983;
t1147 = m(1) + m(2) + m(3);
t975 = t976 ^ 2;
t1447 = -t1035 * t899 - t1147 * t907 - t895 * t1337 + ((-t1140 * t987 - t1249 * (t987 + t975)) * t1158 + (t1124 * t988 + t1233 * t1451) * t1382) * t1154 - t975 * t1156 * t1233;
t978 = t981 ^ 2;
t1446 = -t1044 * t900 - t1147 * t908 - t896 * t1336 + ((-t1140 * t989 - t1248 * (t989 + t978)) * t1167 + (t1124 * t992 + t1232 * t1450) * t1381) * t1154 - t978 * t1156 * t1232;
t979 = t982 ^ 2;
t1445 = -t1045 * t901 - t1147 * t909 - t897 * t1335 + ((-t1140 * t990 - t1247 * (t990 + t979)) * t1169 + (t1124 * t993 + t1231 * t1449) * t1380) * t1154 - t979 * t1156 * t1231;
t980 = t983 ^ 2;
t1444 = -t1046 * t902 - t1147 * t910 - t898 * t1334 + ((-t1140 * t991 - t1246 * (t991 + t980)) * t1171 + (t1124 * t994 + t1230 * t1448) * t1379) * t1154 - t980 * t1156 * t1230;
t1439 = t1023 * t1447;
t1438 = t1024 * t1446;
t1437 = t1025 * t1445;
t1436 = t1026 * t1444;
t1435 = t1465 * t1022 + t1458 * t1040;
t1434 = t1464 * t1021 + t1459 * t1039;
t1433 = -t1463 * t1020 + t1461 * t1038;
t1432 = -t1462 * t1016 + t1460 * t1036;
t1152 = t1180 ^ 2;
t1420 = pkin(2) * t1157;
t1419 = pkin(2) * t1166;
t1418 = pkin(2) * t1168;
t1417 = pkin(2) * t1170;
t1400 = Ifges(3,3) * t1199;
t1361 = t1016 * t1199;
t1360 = t1020 * t1199;
t1359 = t1021 * t1199;
t1358 = t1022 * t1199;
t1357 = t1023 * t1152;
t1355 = t1024 * t1152;
t1353 = t1025 * t1152;
t1351 = t1026 * t1152;
t1349 = t1035 * t1199;
t1348 = t1036 * t1132;
t1347 = t1036 * t1136;
t1346 = t1038 * t1133;
t1345 = t1038 * t1137;
t1344 = t1039 * t1134;
t1343 = t1039 * t1138;
t1342 = t1040 * t1135;
t1341 = t1040 * t1139;
t1340 = t1044 * t1199;
t1339 = t1045 * t1199;
t1338 = t1046 * t1199;
t1333 = t1073 * t1199;
t1332 = t1081 * t1199;
t1331 = t1082 * t1199;
t1330 = t1083 * t1199;
t1328 = t1154 * t1155;
t1257 = t1132 * t1361;
t1256 = t1136 * t1361;
t1255 = t1133 * t1360;
t1254 = t1137 * t1360;
t1253 = t1134 * t1359;
t1252 = t1138 * t1359;
t1251 = t1135 * t1358;
t1250 = t1139 * t1358;
t1221 = pkin(3) * t1327 - t1095 * t1156;
t1220 = pkin(3) * t1323 - t1097 * t1156;
t1219 = pkin(3) * t1321 - t1098 * t1156;
t1218 = pkin(3) * t1319 - t1099 * t1156;
t1217 = Ifges(3,3) * t1361 + t1036 * t1073;
t1216 = Ifges(3,3) * t1360 + t1038 * t1081;
t1215 = Ifges(3,3) * t1359 + t1039 * t1082;
t1214 = Ifges(3,3) * t1358 + t1040 * t1083;
t1213 = t1016 * t1333 + t1031 * t1036;
t1212 = t1020 * t1332 + t1032 * t1038;
t1211 = t1021 * t1331 + t1033 * t1039;
t1210 = t1022 * t1330 + t1034 * t1040;
t1209 = t1023 * (t1087 * t1136 + t1091 * t1132);
t1208 = t1024 * (t1088 * t1137 + t1092 * t1133);
t1207 = t1025 * (t1089 * t1138 + t1093 * t1134);
t1206 = t1026 * (t1090 * t1139 + t1094 * t1135);
t1204 = t1016 * t1349 + t1036 * t1337;
t1203 = t1020 * t1340 + t1038 * t1336;
t1202 = t1021 * t1339 + t1039 * t1335;
t1201 = t1022 * t1338 + t1040 * t1334;
t1188 = mrSges(4,1);
t1187 = mrSges(4,2);
t1030 = -t1102 * t1153 + t1218 * t1155;
t1029 = -t1101 * t1153 + t1219 * t1155;
t1028 = -t1100 * t1153 + t1220 * t1155;
t1027 = -t1096 * t1153 + t1221 * t1155;
t1014 = -t1076 * t1413 + t1102 * t1307 + (pkin(2) * t1319 + t1218 * t1176) * t1153;
t1013 = -t1075 * t1414 + t1101 * t1308 + (pkin(2) * t1321 + t1219 * t1174) * t1153;
t1012 = -t1074 * t1415 + t1100 * t1309 + (pkin(2) * t1323 + t1220 * t1172) * t1153;
t1011 = -t1071 * t1416 + t1096 * t1310 + (pkin(2) * t1327 + t1221 * t1159) * t1153;
t1010 = -(t1079 * t1135 - t1139 * t1318) * t1413 + (t1030 * t1135 + t1058 * t1139) * t1176 + (t1135 * t1328 + t1139 * t1156) * t1417;
t1009 = -(t1078 * t1134 - t1138 * t1320) * t1414 + (t1029 * t1134 + t1057 * t1138) * t1174 + (t1134 * t1328 + t1138 * t1156) * t1418;
t1008 = -(t1077 * t1133 - t1137 * t1322) * t1415 + (t1028 * t1133 + t1056 * t1137) * t1172 + (t1133 * t1328 + t1137 * t1156) * t1419;
t1007 = (t1079 * t1139 + t1135 * t1318) * t1413 + (-t1030 * t1139 + t1058 * t1135) * t1176 + (t1135 * t1156 - t1139 * t1328) * t1417;
t1006 = (t1078 * t1138 + t1134 * t1320) * t1414 + (-t1029 * t1138 + t1057 * t1134) * t1174 + (t1134 * t1156 - t1138 * t1328) * t1418;
t1005 = (t1077 * t1137 + t1133 * t1322) * t1415 + (-t1028 * t1137 + t1056 * t1133) * t1172 + (t1133 * t1156 - t1137 * t1328) * t1419;
t1004 = -(t1072 * t1132 - t1136 * t1326) * t1416 + (t1027 * t1132 + t1055 * t1136) * t1159 + (t1132 * t1328 + t1136 * t1156) * t1420;
t1003 = (t1072 * t1136 + t1132 * t1326) * t1416 + (-t1027 * t1136 + t1055 * t1132) * t1159 + (t1132 * t1156 - t1136 * t1328) * t1420;
t1002 = t1040 * t1206;
t1001 = t1039 * t1207;
t1000 = t1038 * t1208;
t999 = t1036 * t1209;
t998 = t1206 * t1358;
t997 = t1207 * t1359;
t996 = t1208 * t1360;
t995 = t1209 * t1361;
t986 = (t1014 * t1046 + t1019 * t1400 + t1043 * t1083) * t1026;
t985 = (t1013 * t1045 + t1018 * t1400 + t1042 * t1082) * t1025;
t984 = (t1012 * t1044 + t1017 * t1400 + t1041 * t1081) * t1024;
t977 = (t1011 * t1035 + t1015 * t1400 + t1037 * t1073) * t1023;
t974 = (t1014 * t1147 + t1019 * t1338 + t1043 * t1334) * t1026;
t973 = (t1013 * t1147 + t1018 * t1339 + t1042 * t1335) * t1025;
t972 = (t1012 * t1147 + t1017 * t1340 + t1041 * t1336) * t1024;
t971 = (t1011 * t1147 + t1015 * t1349 + t1037 * t1337) * t1023;
t970 = (t1014 * t1334 + t1019 * t1330 + t1034 * t1043) * t1026;
t969 = (t1013 * t1335 + t1018 * t1331 + t1033 * t1042) * t1025;
t968 = (t1012 * t1336 + t1017 * t1332 + t1032 * t1041) * t1024;
t967 = (t1011 * t1337 + t1015 * t1333 + t1031 * t1037) * t1023;
t966 = (-t1007 * t1090 + t1010 * t1094) * t1026;
t965 = (-t1006 * t1089 + t1009 * t1093) * t1025;
t964 = (-t1005 * t1088 + t1008 * t1092) * t1024;
t963 = (-t1003 * t1087 + t1004 * t1091) * t1023;
t962 = (t1010 * t1046 + t1214 * t1135) * t1026;
t961 = (t1009 * t1045 + t1215 * t1134) * t1025;
t960 = (t1008 * t1044 + t1216 * t1133) * t1024;
t959 = (t1007 * t1046 - t1214 * t1139) * t1026;
t958 = (t1006 * t1045 - t1215 * t1138) * t1025;
t957 = (t1005 * t1044 - t1216 * t1137) * t1024;
t956 = (t1004 * t1035 + t1217 * t1132) * t1023;
t955 = (t1003 * t1035 - t1217 * t1136) * t1023;
t954 = (t1010 * t1147 + t1201 * t1135) * t1026;
t953 = (t1009 * t1147 + t1202 * t1134) * t1025;
t952 = (t1008 * t1147 + t1203 * t1133) * t1024;
t951 = (t1007 * t1147 - t1201 * t1139) * t1026;
t950 = (t1006 * t1147 - t1202 * t1138) * t1025;
t949 = (t1005 * t1147 - t1203 * t1137) * t1024;
t948 = (t1004 * t1147 + t1204 * t1132) * t1023;
t947 = (t1003 * t1147 - t1204 * t1136) * t1023;
t946 = (t1010 * t1334 + t1210 * t1135) * t1026;
t945 = (t1009 * t1335 + t1211 * t1134) * t1025;
t944 = (t1008 * t1336 + t1212 * t1133) * t1024;
t943 = (t1007 * t1334 - t1210 * t1139) * t1026;
t942 = (t1006 * t1335 - t1211 * t1138) * t1025;
t941 = (t1005 * t1336 - t1212 * t1137) * t1024;
t940 = (t1004 * t1337 + t1213 * t1132) * t1023;
t939 = (t1003 * t1337 - t1213 * t1136) * t1023;
t930 = Ifges(3,3) * t998 + t1002 * t1083 + t1046 * t966;
t929 = Ifges(3,3) * t997 + t1001 * t1082 + t1045 * t965;
t928 = Ifges(3,3) * t996 + t1000 * t1081 + t1044 * t964;
t927 = t1002 * t1334 + t1046 * t998 + t1147 * t966;
t926 = t1001 * t1335 + t1045 * t997 + t1147 * t965;
t925 = t1000 * t1336 + t1044 * t996 + t1147 * t964;
t924 = Ifges(3,3) * t995 + t1035 * t963 + t1073 * t999;
t923 = t1035 * t995 + t1147 * t963 + t999 * t1337;
t922 = t1002 * t1034 + t1083 * t998 + t966 * t1334;
t921 = t1001 * t1033 + t1082 * t997 + t965 * t1335;
t920 = t1000 * t1032 + t1081 * t996 + t964 * t1336;
t919 = t1031 * t999 + t1073 * t995 + t963 * t1337;
t1 = [(-(t1003 * t947 - t955 * t1256 - t939 * t1347) * t1091 - (t1004 * t947 + t955 * t1257 + t939 * t1348) * t1087) * t1357 + (-(t1007 * t951 - t959 * t1250 - t943 * t1341) * t1094 - (t1010 * t951 + t959 * t1251 + t943 * t1342) * t1090) * t1351 + (-(t1006 * t950 - t958 * t1252 - t942 * t1343) * t1093 - (t1009 * t950 + t958 * t1253 + t942 * t1344) * t1089) * t1353 + (-(t1005 * t949 - t957 * t1254 - t941 * t1345) * t1092 - (t1008 * t949 + t957 * t1255 + t941 * t1346) * t1088) * t1355 + t1152 * (t1144 * t1187 - t1145 * t1188) + t1003 * t1439 + t1005 * t1438 + t1006 * t1437 + t1007 * t1436 - t1435 * t1139 - t1434 * t1138 + t1433 * t1137 + t1432 * t1136; (-(t1003 * t948 - t1256 * t956 - t1347 * t940) * t1091 - (t1004 * t948 + t1257 * t956 + t1348 * t940) * t1087) * t1357 + (-(t1007 * t954 - t1250 * t962 - t1341 * t946) * t1094 - (t1010 * t954 + t1251 * t962 + t1342 * t946) * t1090) * t1351 + (-(t1006 * t953 - t1252 * t961 - t1343 * t945) * t1093 - (t1009 * t953 + t1253 * t961 + t1344 * t945) * t1089) * t1353 + (-(t1005 * t952 - t1254 * t960 - t1345 * t944) * t1092 - (t1008 * t952 + t1255 * t960 + t1346 * t944) * t1088) * t1355 - t1152 * (t1144 * t1188 + t1145 * t1187) + t1004 * t1439 + t1008 * t1438 + t1009 * t1437 + t1010 * t1436 + t1435 * t1135 + t1434 * t1134 - t1433 * t1133 - t1432 * t1132; (-(t1007 * t974 - t1250 * t986 - t1341 * t970) * t1094 - (t1010 * t974 + t1251 * t986 + t1342 * t970) * t1090) * t1351 + (-(t1006 * t973 - t1252 * t985 - t1343 * t969) * t1093 - (t1009 * t973 + t1253 * t985 + t1344 * t969) * t1089) * t1353 + (-(t1005 * t972 - t1254 * t984 - t1345 * t968) * t1092 - (t1008 * t972 + t1255 * t984 + t1346 * t968) * t1088) * t1355 + (-(t1003 * t971 - t1256 * t977 - t1347 * t967) * t1091 - (t1004 * t971 + t1257 * t977 + t1348 * t967) * t1087) * t1357 + t1011 * t1439 + t1012 * t1438 + t1013 * t1437 + t1014 * t1436 + t1458 * t1043 + t1459 * t1042 - t1461 * t1041 - t1460 * t1037 + t1465 * t1019 + t1464 * t1018 + t1463 * t1017 + t1462 * t1015; (-(t1003 * t923 - t1256 * t924 - t1347 * t919) * t1091 - (t1004 * t923 + t1257 * t924 + t1348 * t919) * t1087) * t1357 + (-(t1007 * t927 - t1250 * t930 - t1341 * t922) * t1094 - (t1010 * t927 + t1251 * t930 + t1342 * t922) * t1090) * t1351 + (-(t1006 * t926 - t1252 * t929 - t1343 * t921) * t1093 - (t1009 * t926 + t1253 * t929 + t1344 * t921) * t1089) * t1353 + (-(t1005 * t925 - t1254 * t928 - t1345 * t920) * t1092 - (t1008 * t925 + t1255 * t928 + t1346 * t920) * t1088) * t1355 + (t911 * t1451 + t883) * t999 + t1469 * t998 + t1468 * t997 + t1466 * t996 + t1467 * t995 + t1444 * t966 + t1445 * t965 + t1446 * t964 + t1447 * t963 + (t914 * t1448 + t888) * t1002 + (t913 * t1449 + t887) * t1001 + (t912 * t1450 + t886) * t1000;];
taucX  = t1;
