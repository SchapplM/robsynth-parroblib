% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RRPRR12V1G2A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x15]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:07
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR12V1G2A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:06:14
% EndTime: 2020-08-06 19:06:24
% DurationCPUTime: 10.38s
% Computational Cost: add. (56922->537), mult. (98895->1065), div. (5508->12), fcn. (61494->18), ass. (0->473)
t1205 = legFrame(3,2);
t1192 = sin(t1205);
t1567 = -0.2e1 * t1192;
t1206 = legFrame(2,2);
t1193 = sin(t1206);
t1566 = -0.2e1 * t1193;
t1207 = legFrame(1,2);
t1194 = sin(t1207);
t1565 = -0.2e1 * t1194;
t1195 = cos(t1205);
t1564 = 0.2e1 * t1195;
t1196 = cos(t1206);
t1563 = 0.2e1 * t1196;
t1197 = cos(t1207);
t1562 = 0.2e1 * t1197;
t1221 = xDP(2);
t1223 = pkin(1) + pkin(2);
t1184 = t1221 * t1223;
t1222 = xDP(1);
t1185 = t1222 * t1223;
t1213 = sin(qJ(1,1));
t1219 = cos(qJ(1,1));
t1220 = xDP(3);
t1411 = t1220 * t1223;
t1527 = qJ(3,1) * t1222;
t1528 = qJ(3,1) * t1221;
t1117 = (t1213 * t1185 - t1528) * t1197 + (-t1213 * t1184 - t1527) * t1194 + t1219 * t1411;
t1212 = sin(qJ(2,1));
t1416 = t1212 * t1219;
t1532 = t1213 * pkin(4);
t1329 = qJ(3,1) * t1416 - t1532;
t1150 = t1220 * t1329;
t1529 = t1219 * pkin(4);
t1183 = t1222 * t1529;
t1218 = cos(qJ(2,1));
t1203 = t1218 ^ 2;
t1384 = t1221 * t1529;
t1550 = -t1213 * t1528 + t1185;
t1553 = t1213 * t1527 + t1184;
t1099 = t1117 * t1203 + ((t1553 * t1212 + t1183) * t1197 + (t1550 * t1212 - t1384) * t1194 + t1150) * t1218 + qJ(3,1) * (t1194 * t1222 + t1197 * t1221);
t1191 = t1212 * qJ(3,1);
t1405 = t1223 * t1218;
t1171 = t1191 + t1405;
t1163 = 0.1e1 / t1171;
t1495 = t1099 * t1163;
t1561 = 0.2e1 * t1495;
t1211 = sin(qJ(1,2));
t1217 = cos(qJ(1,2));
t1525 = qJ(3,2) * t1222;
t1526 = qJ(3,2) * t1221;
t1116 = (t1211 * t1185 - t1526) * t1196 + (-t1211 * t1184 - t1525) * t1193 + t1217 * t1411;
t1210 = sin(qJ(2,2));
t1421 = t1210 * t1217;
t1533 = t1211 * pkin(4);
t1330 = qJ(3,2) * t1421 - t1533;
t1149 = t1220 * t1330;
t1530 = t1217 * pkin(4);
t1182 = t1222 * t1530;
t1216 = cos(qJ(2,2));
t1202 = t1216 ^ 2;
t1385 = t1221 * t1530;
t1551 = -t1211 * t1526 + t1185;
t1554 = t1211 * t1525 + t1184;
t1098 = t1116 * t1202 + ((t1554 * t1210 + t1182) * t1196 + (t1551 * t1210 - t1385) * t1193 + t1149) * t1216 + qJ(3,2) * (t1193 * t1222 + t1196 * t1221);
t1190 = t1210 * qJ(3,2);
t1406 = t1223 * t1216;
t1170 = t1190 + t1406;
t1160 = 0.1e1 / t1170;
t1498 = t1098 * t1160;
t1560 = 0.2e1 * t1498;
t1209 = sin(qJ(1,3));
t1215 = cos(qJ(1,3));
t1523 = qJ(3,3) * t1222;
t1524 = qJ(3,3) * t1221;
t1115 = (t1209 * t1185 - t1524) * t1195 + (-t1209 * t1184 - t1523) * t1192 + t1215 * t1411;
t1208 = sin(qJ(2,3));
t1426 = t1208 * t1215;
t1534 = t1209 * pkin(4);
t1331 = qJ(3,3) * t1426 - t1534;
t1148 = t1220 * t1331;
t1531 = t1215 * pkin(4);
t1181 = t1222 * t1531;
t1214 = cos(qJ(2,3));
t1201 = t1214 ^ 2;
t1386 = t1221 * t1531;
t1552 = -t1209 * t1524 + t1185;
t1555 = t1209 * t1523 + t1184;
t1097 = t1115 * t1201 + ((t1555 * t1208 + t1181) * t1195 + (t1552 * t1208 - t1386) * t1192 + t1148) * t1214 + qJ(3,3) * (t1192 * t1222 + t1195 * t1221);
t1189 = t1208 * qJ(3,3);
t1407 = t1223 * t1214;
t1169 = t1189 + t1407;
t1157 = 0.1e1 / t1169;
t1501 = t1097 * t1157;
t1559 = 0.2e1 * t1501;
t1105 = -t1194 * t1384 + t1117 * t1218 + t1183 * t1197 + t1150 + (t1550 * t1194 + t1553 * t1197) * t1212;
t1231 = 0.1e1 / qJ(3,1);
t1102 = t1105 * t1231;
t1408 = t1223 * t1212;
t1514 = t1218 * qJ(3,1);
t1280 = -t1408 + t1514;
t1135 = (-t1194 * t1221 + t1197 * t1222) * t1219 - t1213 * t1220;
t1164 = 0.1e1 / t1171 ^ 2;
t1165 = t1163 * t1164;
t1439 = t1165 * t1231;
t1299 = t1099 * t1135 * t1439;
t1469 = t1135 * t1231;
t1351 = t1164 * t1469;
t1402 = t1223 * t1231;
t1471 = t1135 * t1164;
t1538 = pkin(4) * t1135;
t1045 = -(t1102 * t1212 + (-t1538 + (-t1212 * t1402 + t1218) * t1099) * t1163) * t1471 - t1280 * t1299 - t1212 * t1105 * t1351;
t1188 = 0.2e1 * t1203 - 0.1e1;
t1412 = t1218 * t1231;
t1337 = t1212 * t1412;
t1549 = -2 * pkin(1);
t1054 = t1105 * t1337 + (t1337 * t1549 + t1188) * t1495;
t1177 = t1218 * pkin(1) + t1191;
t1441 = t1163 * t1218;
t1558 = 0.2e1 * (t1045 * t1177 * t1441 + t1054 * t1471) * t1219;
t1104 = -t1193 * t1385 + t1116 * t1216 + t1182 * t1196 + t1149 + (t1551 * t1193 + t1554 * t1196) * t1210;
t1228 = 0.1e1 / qJ(3,2);
t1101 = t1104 * t1228;
t1409 = t1223 * t1210;
t1515 = t1216 * qJ(3,2);
t1279 = -t1409 + t1515;
t1134 = (-t1193 * t1221 + t1196 * t1222) * t1217 - t1211 * t1220;
t1161 = 0.1e1 / t1170 ^ 2;
t1162 = t1160 * t1161;
t1445 = t1162 * t1228;
t1300 = t1098 * t1134 * t1445;
t1473 = t1134 * t1228;
t1353 = t1161 * t1473;
t1403 = t1223 * t1228;
t1475 = t1134 * t1161;
t1539 = pkin(4) * t1134;
t1044 = -(t1101 * t1210 + (-t1539 + (-t1210 * t1403 + t1216) * t1098) * t1160) * t1475 - t1279 * t1300 - t1210 * t1104 * t1353;
t1187 = 0.2e1 * t1202 - 0.1e1;
t1413 = t1216 * t1228;
t1339 = t1210 * t1413;
t1053 = t1104 * t1339 + (t1339 * t1549 + t1187) * t1498;
t1175 = t1216 * pkin(1) + t1190;
t1447 = t1160 * t1216;
t1557 = 0.2e1 * (t1044 * t1175 * t1447 + t1053 * t1475) * t1217;
t1103 = -t1192 * t1386 + t1115 * t1214 + t1181 * t1195 + t1148 + (t1552 * t1192 + t1555 * t1195) * t1208;
t1225 = 0.1e1 / qJ(3,3);
t1100 = t1103 * t1225;
t1410 = t1223 * t1208;
t1516 = t1214 * qJ(3,3);
t1278 = -t1410 + t1516;
t1133 = (-t1192 * t1221 + t1195 * t1222) * t1215 - t1209 * t1220;
t1158 = 0.1e1 / t1169 ^ 2;
t1159 = t1157 * t1158;
t1451 = t1159 * t1225;
t1301 = t1097 * t1133 * t1451;
t1477 = t1133 * t1225;
t1355 = t1158 * t1477;
t1404 = t1223 * t1225;
t1479 = t1133 * t1158;
t1540 = pkin(4) * t1133;
t1043 = -(t1100 * t1208 + (-t1540 + (-t1208 * t1404 + t1214) * t1097) * t1157) * t1479 - t1278 * t1301 - t1208 * t1103 * t1355;
t1186 = 0.2e1 * t1201 - 0.1e1;
t1414 = t1214 * t1225;
t1341 = t1208 * t1414;
t1052 = t1103 * t1341 + (t1341 * t1549 + t1186) * t1501;
t1173 = t1214 * pkin(1) + t1189;
t1453 = t1157 * t1214;
t1556 = 0.2e1 * (t1043 * t1173 * t1453 + t1052 * t1479) * t1215;
t1224 = qJ(3,3) ^ 2;
t1233 = pkin(4) ^ 2;
t1480 = t1133 * t1157;
t1356 = t1201 * t1480;
t1452 = t1157 * t1225;
t1361 = t1097 * t1452;
t1089 = pkin(1) * t1361;
t1070 = t1089 - t1100;
t1061 = pkin(2) * t1361 + t1070;
t1429 = t1208 * t1061;
t1500 = t1097 * t1225;
t1541 = -pkin(4) / 0.2e1;
t1046 = (qJ(3,3) + t1223) * (-qJ(3,3) + t1223) * t1356 + 0.2e1 * (t1133 * t1410 + t1500 * t1541) * qJ(3,3) * t1453 + pkin(4) * t1429 + (t1224 + t1233) * t1480;
t1455 = t1157 * t1208;
t1124 = t1455 * t1540;
t1082 = t1223 * t1361 + t1124;
t1235 = pkin(1) ^ 2;
t1335 = -t1235 + (t1549 - pkin(2)) * pkin(2);
t1226 = 0.1e1 / qJ(3,3) ^ 2;
t1492 = t1103 * t1226;
t1499 = t1097 * t1226;
t1244 = (t1046 * t1477 - (t1103 * t1404 + ((-t1224 + t1335) * t1500 + t1278 * t1540) * t1157) * t1499) * t1157 - t1082 * t1492;
t1227 = qJ(3,2) ^ 2;
t1476 = t1134 * t1160;
t1354 = t1202 * t1476;
t1446 = t1160 * t1228;
t1359 = t1098 * t1446;
t1091 = pkin(1) * t1359;
t1071 = t1091 - t1101;
t1062 = pkin(2) * t1359 + t1071;
t1424 = t1210 * t1062;
t1497 = t1098 * t1228;
t1047 = (qJ(3,2) + t1223) * (-qJ(3,2) + t1223) * t1354 + 0.2e1 * qJ(3,2) * (t1134 * t1409 + t1497 * t1541) * t1447 + pkin(4) * t1424 + (t1227 + t1233) * t1476;
t1449 = t1160 * t1210;
t1125 = t1449 * t1539;
t1083 = t1223 * t1359 + t1125;
t1229 = 0.1e1 / qJ(3,2) ^ 2;
t1491 = t1104 * t1229;
t1496 = t1098 * t1229;
t1243 = (t1047 * t1473 - (t1104 * t1403 + ((-t1227 + t1335) * t1497 + t1279 * t1539) * t1160) * t1496) * t1160 - t1083 * t1491;
t1230 = qJ(3,1) ^ 2;
t1472 = t1135 * t1163;
t1352 = t1203 * t1472;
t1440 = t1163 * t1231;
t1357 = t1099 * t1440;
t1093 = pkin(1) * t1357;
t1072 = t1093 - t1102;
t1063 = pkin(2) * t1357 + t1072;
t1419 = t1212 * t1063;
t1494 = t1099 * t1231;
t1048 = (qJ(3,1) + t1223) * (-qJ(3,1) + t1223) * t1352 + 0.2e1 * qJ(3,1) * (t1135 * t1408 + t1494 * t1541) * t1441 + pkin(4) * t1419 + (t1230 + t1233) * t1472;
t1443 = t1163 * t1212;
t1126 = t1443 * t1538;
t1084 = t1223 * t1357 + t1126;
t1232 = 0.1e1 / qJ(3,1) ^ 2;
t1490 = t1105 * t1232;
t1493 = t1099 * t1232;
t1242 = (t1048 * t1469 - (t1105 * t1402 + ((-t1230 + t1335) * t1494 + t1280 * t1538) * t1163) * t1493) * t1163 - t1084 * t1490;
t1548 = 2 * pkin(1);
t1547 = -0.2e1 * t1201;
t1546 = -0.2e1 * t1202;
t1545 = -0.2e1 * t1203;
t1544 = -0.2e1 * t1209;
t1543 = -0.2e1 * t1211;
t1542 = -0.2e1 * t1213;
t1085 = t1097 * t1455;
t1037 = t1214 * t1046 * t1355 - (-(t1124 + t1061) * t1407 + (pkin(4) * t1356 - t1429) * qJ(3,3)) * t1158 * t1499 - (t1082 * t1214 + t1085) * t1157 * t1492;
t1537 = t1037 * pkin(1);
t1086 = t1098 * t1449;
t1038 = t1216 * t1047 * t1353 - (-(t1125 + t1062) * t1406 + (pkin(4) * t1354 - t1424) * qJ(3,2)) * t1161 * t1496 - (t1083 * t1216 + t1086) * t1160 * t1491;
t1536 = t1038 * pkin(1);
t1087 = t1099 * t1443;
t1039 = t1218 * t1048 * t1351 - (-(t1126 + t1063) * t1405 + (pkin(4) * t1352 - t1419) * qJ(3,1)) * t1164 * t1493 - (t1084 * t1218 + t1087) * t1163 * t1490;
t1535 = t1039 * pkin(1);
t1522 = t1192 * qJ(3,3);
t1521 = t1193 * qJ(3,2);
t1520 = t1194 * qJ(3,1);
t1519 = t1195 * qJ(3,3);
t1518 = t1196 * qJ(3,2);
t1517 = t1197 * qJ(3,1);
t1513 = (-t1244 + 0.2e1 * t1537) * t1157;
t1512 = (-t1243 + 0.2e1 * t1536) * t1160;
t1511 = (-t1242 + 0.2e1 * t1535) * t1163;
t1510 = t1037 * t1214;
t1509 = t1038 * t1216;
t1508 = t1039 * t1218;
t1130 = t1133 ^ 2;
t1127 = t1130 * t1158;
t1483 = t1130 * t1201;
t1504 = t1097 ^ 2 * t1226;
t1076 = -t1127 + (t1483 - t1504) * t1158;
t1507 = t1076 * t1225;
t1131 = t1134 ^ 2;
t1128 = t1131 * t1161;
t1482 = t1131 * t1202;
t1503 = t1098 ^ 2 * t1229;
t1077 = -t1128 + (t1482 - t1503) * t1161;
t1506 = t1077 * t1228;
t1132 = t1135 ^ 2;
t1129 = t1132 * t1164;
t1481 = t1132 * t1203;
t1502 = t1099 ^ 2 * t1232;
t1078 = -t1129 + (t1481 - t1502) * t1164;
t1505 = t1078 * t1231;
t1428 = t1208 * t1209;
t1154 = qJ(3,3) * t1428 + t1531;
t1425 = t1209 * t1223;
t1109 = (t1195 * t1425 - t1522) * t1201 + (t1154 * t1195 + t1192 * t1410) * t1214 + t1522;
t1489 = t1109 * t1225;
t1423 = t1210 * t1211;
t1155 = qJ(3,2) * t1423 + t1530;
t1420 = t1211 * t1223;
t1110 = (t1196 * t1420 - t1521) * t1202 + (t1155 * t1196 + t1193 * t1409) * t1216 + t1521;
t1488 = t1110 * t1228;
t1418 = t1212 * t1213;
t1156 = qJ(3,1) * t1418 + t1529;
t1415 = t1213 * t1223;
t1111 = (t1197 * t1415 - t1520) * t1203 + (t1156 * t1197 + t1194 * t1408) * t1218 + t1520;
t1487 = t1111 * t1231;
t1112 = (-t1192 * t1425 - t1519) * t1201 + (-t1192 * t1154 + t1195 * t1410) * t1214 + t1519;
t1486 = t1112 * t1225;
t1113 = (-t1193 * t1420 - t1518) * t1202 + (-t1193 * t1155 + t1196 * t1409) * t1216 + t1518;
t1485 = t1113 * t1228;
t1114 = (-t1194 * t1415 - t1517) * t1203 + (-t1194 * t1156 + t1197 * t1408) * t1218 + t1517;
t1484 = t1114 * t1231;
t1478 = t1133 * t1215;
t1474 = t1134 * t1217;
t1470 = t1135 * t1219;
t1139 = t1215 * t1407 + t1331;
t1468 = t1139 * t1214;
t1140 = t1217 * t1406 + t1330;
t1467 = t1140 * t1216;
t1141 = t1219 * t1405 + t1329;
t1466 = t1141 * t1218;
t1142 = t1169 * t1215 - t1534;
t1465 = t1142 * t1208;
t1464 = t1142 * t1225;
t1143 = t1170 * t1217 - t1533;
t1463 = t1143 * t1210;
t1462 = t1143 * t1228;
t1144 = t1171 * t1219 - t1532;
t1461 = t1144 * t1212;
t1460 = t1144 * t1231;
t1427 = t1208 * t1214;
t1334 = t1427 * t1548;
t1392 = qJ(3,3) * t1547;
t1459 = (qJ(3,3) + t1334 + t1392) * t1159;
t1422 = t1210 * t1216;
t1333 = t1422 * t1548;
t1391 = qJ(3,2) * t1546;
t1458 = (qJ(3,2) + t1333 + t1391) * t1162;
t1417 = t1212 * t1218;
t1332 = t1417 * t1548;
t1390 = qJ(3,1) * t1545;
t1457 = (qJ(3,1) + t1332 + t1390) * t1165;
t1172 = -t1208 * pkin(1) + t1516;
t1456 = t1157 * t1172;
t1454 = t1157 * t1209;
t1174 = -t1210 * pkin(1) + t1515;
t1450 = t1160 * t1174;
t1448 = t1160 * t1211;
t1176 = -t1212 * pkin(1) + t1514;
t1444 = t1163 * t1176;
t1442 = t1163 * t1213;
t1438 = t1192 * t1215;
t1437 = t1193 * t1217;
t1436 = t1194 * t1219;
t1435 = t1195 * t1215;
t1434 = t1196 * t1217;
t1433 = t1197 * t1219;
t1432 = t1201 * t1225;
t1431 = t1202 * t1228;
t1430 = t1203 * t1231;
t1032 = t1244 - t1537;
t1298 = t1427 * t1127;
t1401 = -pkin(1) * t1298 + t1076 * qJ(3,3) + t1032;
t1034 = t1243 - t1536;
t1297 = t1422 * t1128;
t1400 = -pkin(1) * t1297 + t1077 * qJ(3,2) + t1034;
t1036 = t1242 - t1535;
t1296 = t1417 * t1129;
t1399 = -pkin(1) * t1296 + t1078 * qJ(3,1) + t1036;
t1398 = -t1224 + t1235;
t1397 = -t1227 + t1235;
t1396 = -t1230 + t1235;
t1389 = t1133 * t1544;
t1388 = t1134 * t1543;
t1387 = t1135 * t1542;
t1377 = t1037 * t1452;
t1376 = t1038 * t1446;
t1375 = t1039 * t1440;
t1374 = t1043 * t1426;
t1373 = t1044 * t1421;
t1372 = t1045 * t1416;
t1371 = (t1085 + (t1089 - 0.2e1 * t1100) * t1214) * t1097 * t1158;
t1370 = (t1086 + (t1091 - 0.2e1 * t1101) * t1216) * t1098 * t1161;
t1369 = (t1087 + (t1093 - 0.2e1 * t1102) * t1218) * t1099 * t1164;
t1365 = t1159 * t1504;
t1364 = t1162 * t1503;
t1363 = t1165 * t1502;
t1362 = t1097 * t1478;
t1360 = t1098 * t1474;
t1358 = t1099 * t1470;
t1350 = t1139 * t1414;
t1349 = t1140 * t1413;
t1348 = t1141 * t1412;
t1347 = t1157 * t1438;
t1346 = t1157 * t1435;
t1345 = t1160 * t1437;
t1344 = t1160 * t1434;
t1343 = t1163 * t1436;
t1342 = t1163 * t1433;
t1340 = t1043 * t1454;
t1338 = t1044 * t1448;
t1336 = t1045 * t1442;
t1322 = t1103 * t1559;
t1328 = ((t1224 + t1235) * t1037 - pkin(1) * t1244 + t1225 * t1322 + (-t1186 * qJ(3,3) * pkin(1) + t1398 * t1427) * t1127) * t1157;
t1321 = t1104 * t1560;
t1327 = ((t1227 + t1235) * t1038 - pkin(1) * t1243 + t1228 * t1321 + (-t1187 * qJ(3,2) * pkin(1) + t1397 * t1422) * t1128) * t1160;
t1320 = t1105 * t1561;
t1326 = ((t1230 + t1235) * t1039 - pkin(1) * t1242 + t1231 * t1320 + (-t1188 * qJ(3,1) * pkin(1) + t1396 * t1417) * t1129) * t1163;
t1325 = -0.2e1 * t1043 * t1428;
t1324 = -0.2e1 * t1044 * t1423;
t1323 = -0.2e1 * t1045 * t1418;
t1319 = t1158 * t1389;
t1318 = t1161 * t1388;
t1317 = t1164 * t1387;
t1316 = (t1226 * t1322 + (pkin(1) * t1547 - 0.2e1 * qJ(3,3) * t1427 + pkin(1)) * t1127) * t1225 + 0.2e1 * t1037;
t1315 = (t1229 * t1321 + (pkin(1) * t1546 - 0.2e1 * qJ(3,2) * t1422 + pkin(1)) * t1128) * t1228 + 0.2e1 * t1038;
t1314 = (t1232 * t1320 + (pkin(1) * t1545 - 0.2e1 * qJ(3,1) * t1417 + pkin(1)) * t1129) * t1231 + 0.2e1 * t1039;
t1049 = -t1103 * t1432 + (pkin(1) * t1432 + t1427) * t1559 - t1070;
t1313 = t1049 * t1158 * t1478;
t1050 = -t1104 * t1431 + (pkin(1) * t1431 + t1422) * t1560 - t1071;
t1312 = t1050 * t1161 * t1474;
t1051 = -t1105 * t1430 + (pkin(1) * t1430 + t1417) * t1561 - t1072;
t1311 = t1051 * t1164 * t1470;
t1310 = t1215 * t1371;
t1309 = t1217 * t1370;
t1308 = t1219 * t1369;
t1307 = t1209 * t1365;
t1306 = t1215 * t1365;
t1305 = t1211 * t1364;
t1304 = t1217 * t1364;
t1303 = t1213 * t1363;
t1302 = t1219 * t1363;
t1295 = t1159 * t1341;
t1294 = t1162 * t1339;
t1293 = t1165 * t1337;
t1292 = t1043 * t1347;
t1291 = t1044 * t1345;
t1290 = t1045 * t1343;
t1289 = t1043 * t1346;
t1288 = t1044 * t1344;
t1287 = t1045 * t1342;
t1286 = t1374 * t1567;
t1285 = t1374 * t1564;
t1284 = t1373 * t1566;
t1283 = t1373 * t1563;
t1282 = t1372 * t1565;
t1281 = t1372 * t1562;
t1277 = t1208 * t1306;
t1276 = t1214 * t1306;
t1275 = t1210 * t1304;
t1274 = t1216 * t1304;
t1273 = t1212 * t1302;
t1272 = t1218 * t1302;
t1271 = t1186 * t1301;
t1270 = t1187 * t1300;
t1269 = t1188 * t1299;
t1268 = -t1037 - t1298;
t1267 = -t1038 - t1297;
t1266 = -t1039 - t1296;
t1265 = t1215 * t1271;
t1264 = t1217 * t1270;
t1263 = t1219 * t1269;
t1262 = t1130 * t1459 + t1513;
t1261 = t1131 * t1458 + t1512;
t1260 = t1132 * t1457 + t1511;
t1259 = -t1157 * (qJ(3,3) * t1334 + t1398 * t1201 + t1224) * t1043 + 0.2e1 * ((t1089 - t1100 / 0.2e1) * t1392 + (pkin(1) * t1070 - qJ(3,3) * t1501) * t1427 + qJ(3,3) * t1070) * t1479;
t1258 = -t1160 * (qJ(3,2) * t1333 + t1397 * t1202 + t1227) * t1044 + 0.2e1 * ((t1091 - t1101 / 0.2e1) * t1391 + (pkin(1) * t1071 - qJ(3,2) * t1498) * t1422 + qJ(3,2) * t1071) * t1475;
t1257 = -t1163 * (qJ(3,1) * t1332 + t1396 * t1203 + t1230) * t1045 + 0.2e1 * ((t1093 - t1102 / 0.2e1) * t1390 + (pkin(1) * t1072 - qJ(3,1) * t1495) * t1417 + qJ(3,1) * t1072) * t1471;
t1250 = t1157 * (t1037 * t1435 + t1043 * t1489);
t1249 = t1157 * (-t1037 * t1438 + t1043 * t1486);
t1248 = t1160 * (t1038 * t1434 + t1044 * t1488);
t1247 = t1160 * (-t1038 * t1437 + t1044 * t1485);
t1246 = t1163 * (t1039 * t1433 + t1045 * t1487);
t1245 = t1163 * (-t1039 * t1436 + t1045 * t1484);
t1200 = t1212 ^ 2;
t1199 = t1210 ^ 2;
t1198 = t1208 ^ 2;
t1138 = t1171 * t1213 + t1529;
t1137 = t1170 * t1211 + t1530;
t1136 = t1169 * t1209 + t1531;
t1123 = -t1138 * t1194 - t1197 * t1280;
t1122 = t1138 * t1197 - t1194 * t1280;
t1121 = -t1137 * t1193 - t1196 * t1279;
t1120 = t1137 * t1196 - t1193 * t1279;
t1119 = -t1136 * t1192 - t1195 * t1278;
t1118 = t1136 * t1195 - t1192 * t1278;
t1108 = -0.2e1 * t1164 * t1481 + t1129;
t1107 = -0.2e1 * t1161 * t1482 + t1128;
t1106 = -0.2e1 * t1158 * t1483 + t1127;
t1027 = qJ(3,1) * t1508 + t1036 * t1212;
t1026 = qJ(3,2) * t1509 + t1034 * t1210;
t1025 = qJ(3,3) * t1510 + t1032 * t1208;
t1 = [t1287 + t1288 + t1289, 0, 0, t1198 * t1289 + t1199 * t1288 + t1200 * t1287 + (-t1111 * t1132 + t1358 * t1562) * t1293 + (-t1110 * t1131 + t1360 * t1563) * t1294 + (-t1109 * t1130 + t1362 * t1564) * t1295, t1265 * t1564 + t1264 * t1563 + t1263 * t1562 + (t1108 * t1487 + t1218 * t1281) * t1163 + (t1107 * t1488 + t1216 * t1283) * t1160 + (t1106 * t1489 + t1214 * t1285) * t1157, t1195 * t1276 + t1196 * t1274 + t1197 * t1272 + t1208 * t1250 + t1210 * t1248 + t1212 * t1246, -t1195 * t1277 - t1196 * t1275 - t1197 * t1273 + t1214 * t1250 + t1216 * t1248 + t1218 * t1246, t1109 * t1377 + t1110 * t1376 + t1111 * t1375, 0, 0, t1197 * t1558 + t1196 * t1557 + t1195 * t1556 + (t1111 * t1260 + t1122 * t1266) * t1231 + (t1110 * t1261 + t1120 * t1267) * t1228 + (t1109 * t1262 + t1118 * t1268) * t1225, t1025 * t1346 + t1026 * t1344 + t1027 * t1342 + (-t1197 * t1308 + (t1111 * t1444 + t1122 * t1212) * t1045) * t1231 + (-t1196 * t1309 + (t1110 * t1450 + t1120 * t1210) * t1044) * t1228 + (-t1195 * t1310 + (t1109 * t1456 + t1118 * t1208) * t1043) * t1225, t1313 * t1564 + t1312 * t1563 + t1311 * t1562 + t1118 * t1507 + t1120 * t1506 + t1122 * t1505 + (t1111 * t1314 + t1177 * t1281) * t1163 + (t1110 * t1315 + t1175 * t1283) * t1160 + (t1109 * t1316 + t1173 * t1285) * t1157, -t1257 * t1433 - t1258 * t1434 - t1259 * t1435 + (t1111 * t1326 + t1122 * t1399) * t1231 + (t1110 * t1327 + t1120 * t1400) * t1228 + (t1109 * t1328 + t1118 * t1401) * t1225, 0; -t1290 - t1291 - t1292, 0, 0, -t1198 * t1292 - t1199 * t1291 - t1200 * t1290 + (-t1114 * t1132 + t1358 * t1565) * t1293 + (-t1113 * t1131 + t1360 * t1566) * t1294 + (-t1112 * t1130 + t1362 * t1567) * t1295, t1265 * t1567 + t1264 * t1566 + t1263 * t1565 + (t1108 * t1484 + t1218 * t1282) * t1163 + (t1107 * t1485 + t1216 * t1284) * t1160 + (t1106 * t1486 + t1214 * t1286) * t1157, -t1192 * t1276 - t1193 * t1274 - t1194 * t1272 + t1208 * t1249 + t1210 * t1247 + t1212 * t1245, t1192 * t1277 + t1193 * t1275 + t1194 * t1273 + t1214 * t1249 + t1216 * t1247 + t1218 * t1245, t1112 * t1377 + t1113 * t1376 + t1114 * t1375, 0, 0, -t1194 * t1558 - t1193 * t1557 - t1192 * t1556 + (t1114 * t1260 + t1123 * t1266) * t1231 + (t1113 * t1261 + t1121 * t1267) * t1228 + (t1112 * t1262 + t1119 * t1268) * t1225, -t1025 * t1347 - t1026 * t1345 - t1027 * t1343 + (t1194 * t1308 + (t1114 * t1444 + t1123 * t1212) * t1045) * t1231 + (t1193 * t1309 + (t1113 * t1450 + t1121 * t1210) * t1044) * t1228 + (t1192 * t1310 + (t1112 * t1456 + t1119 * t1208) * t1043) * t1225, t1313 * t1567 + t1312 * t1566 + t1311 * t1565 + t1119 * t1507 + t1121 * t1506 + t1123 * t1505 + (t1114 * t1314 + t1177 * t1282) * t1163 + (t1113 * t1315 + t1175 * t1284) * t1160 + (t1112 * t1316 + t1173 * t1286) * t1157, t1257 * t1436 + t1258 * t1437 + t1259 * t1438 + (t1114 * t1326 + t1123 * t1399) * t1231 + (t1113 * t1327 + t1121 * t1400) * t1228 + (t1112 * t1328 + t1119 * t1401) * t1225, 0; -t1336 - t1338 - t1340, 0, 0, -t1198 * t1340 - t1199 * t1338 - t1200 * t1336 + (t1099 * t1218 * t1387 - t1141 * t1481) * t1212 * t1439 + (t1098 * t1216 * t1388 - t1140 * t1482) * t1210 * t1445 + (t1097 * t1214 * t1389 - t1139 * t1483) * t1208 * t1451, t1271 * t1544 + t1270 * t1543 + t1269 * t1542 + (t1108 * t1141 * t1231 + t1323) * t1441 + (t1107 * t1140 * t1228 + t1324) * t1447 + (t1106 * t1139 * t1225 + t1325) * t1453, -t1214 * t1307 - t1216 * t1305 - t1218 * t1303 + (-t1039 * t1213 + t1045 * t1348) * t1443 + (-t1038 * t1211 + t1044 * t1349) * t1449 + (-t1037 * t1209 + t1043 * t1350) * t1455, t1208 * t1307 + t1210 * t1305 + t1212 * t1303 + (t1141 * t1045 * t1430 - t1213 * t1508) * t1163 + (t1140 * t1044 * t1431 - t1211 * t1509) * t1160 + (t1139 * t1043 * t1432 - t1209 * t1510) * t1157, t1037 * t1157 * t1350 + t1038 * t1160 * t1349 + t1039 * t1163 * t1348, 0, 0, t1052 * t1319 + t1053 * t1318 + t1054 * t1317 - t1037 * t1464 - t1038 * t1462 - t1039 * t1460 + (-0.2e1 * t1177 * t1336 + (t1141 * t1511 + (t1141 * t1457 - t1164 * t1461) * t1132) * t1231) * t1218 + (-0.2e1 * t1175 * t1338 + (t1140 * t1512 + (t1140 * t1458 - t1161 * t1463) * t1131) * t1228) * t1216 + (-0.2e1 * t1173 * t1340 + (t1139 * t1513 + (t1139 * t1459 - t1158 * t1465) * t1130) * t1225) * t1214, -t1025 * t1454 - t1026 * t1448 - t1027 * t1442 + (t1213 * t1369 + (t1141 * t1176 * t1441 + t1461) * t1045) * t1231 + (t1211 * t1370 + (t1140 * t1174 * t1447 + t1463) * t1044) * t1228 + (t1209 * t1371 + (t1139 * t1172 * t1453 + t1465) * t1043) * t1225, t1049 * t1319 + t1050 * t1318 + t1051 * t1317 + t1076 * t1464 + t1077 * t1462 + t1078 * t1460 + (t1177 * t1323 + t1314 * t1466) * t1163 + (t1175 * t1324 + t1315 * t1467) * t1160 + (t1173 * t1325 + t1316 * t1468) * t1157, t1257 * t1213 + t1258 * t1211 + t1259 * t1209 + (t1144 * t1399 + t1326 * t1466) * t1231 + (t1143 * t1400 + t1327 * t1467) * t1228 + (t1142 * t1401 + t1328 * t1468) * t1225, 0;];
tau_reg  = t1;
