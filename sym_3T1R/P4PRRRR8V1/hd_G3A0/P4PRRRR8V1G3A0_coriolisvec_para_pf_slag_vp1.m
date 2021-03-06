% Calculate vector of centrifugal and coriolis load on the joints for
% P4PRRRR8V1G3A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% taucX [4x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-09-20 23:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-09-20 22:58:59
% EndTime: 2020-09-20 22:59:14
% DurationCPUTime: 15.52s
% Computational Cost: add. (38208->480), mult. (109690->981), div. (13200->11), fcn. (114260->38), ass. (0->391)
t1179 = sin(pkin(6));
t1181 = cos(pkin(6));
t1184 = sin(qJ(2,4));
t1186 = cos(qJ(2,4));
t1182 = cos(pkin(3));
t1373 = t1182 * t1186;
t1374 = t1182 * t1184;
t1185 = cos(qJ(3,4));
t1446 = pkin(2) * t1185;
t1069 = (-t1179 * t1184 + t1181 * t1373) * t1446 + pkin(5) * (t1179 * t1186 + t1181 * t1374);
t1070 = (t1179 * t1373 + t1181 * t1184) * t1446 + (t1179 * t1374 - t1181 * t1186) * pkin(5);
t1207 = xDP(3);
t1210 = xP(4);
t1167 = sin(t1210);
t1168 = cos(t1210);
t1217 = koppelP(4,2);
t1221 = koppelP(4,1);
t1125 = t1167 * t1221 + t1168 * t1217;
t1129 = -t1167 * t1217 + t1168 * t1221;
t1187 = legFrame(4,2);
t1159 = sin(t1187);
t1163 = cos(t1187);
t1206 = xDP(4);
t1208 = xDP(2);
t1209 = xDP(1);
t1282 = (t1129 * t1206 + t1208) * t1159 - (-t1125 * t1206 + t1209) * t1163;
t1229 = 0.1e1 / pkin(2);
t1362 = t1184 * t1185;
t1133 = pkin(2) * t1362 - pkin(5) * t1186;
t1180 = sin(pkin(3));
t1183 = sin(qJ(3,4));
t1447 = pkin(2) * t1183;
t1106 = t1133 * t1180 + t1182 * t1447;
t1468 = 0.1e1 / t1106;
t1403 = t1468 / t1185;
t1328 = t1229 * t1403;
t1015 = (t1282 * t1069 + t1070 * t1207) * t1328;
t1246 = t1180 * t1185 + t1183 * t1374;
t1363 = t1183 * t1186;
t1078 = t1246 * t1179 - t1181 * t1363;
t1079 = t1179 * t1363 + t1246 * t1181;
t1031 = (t1078 * t1207 + t1282 * t1079) * t1403;
t1158 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t1213 = 0.2e1 * qJ(3,4);
t1225 = (rSges(3,2) ^ 2);
t1226 = (rSges(3,1) ^ 2);
t1355 = (t1225 + t1226);
t1230 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + ((2 * rSges(3,3) ^ 2 + t1355) * m(3)) / 0.2e1;
t1147 = ((-t1225 + t1226) * m(3)) - Icges(3,1) + Icges(3,2);
t1464 = t1147 / 0.2e1;
t1057 = cos(t1213) * t1464 - t1158 * sin(t1213) + t1230;
t1465 = m(3) * rSges(3,3);
t1156 = rSges(3,2) * t1465 - Icges(3,6);
t1157 = rSges(3,1) * t1465 - Icges(3,5);
t1121 = -t1156 * t1185 - t1157 * t1183;
t1169 = t1185 ^ 2;
t1318 = t1147 * t1183 * t1185;
t1155 = m(2) * rSges(2,2) - t1465;
t1205 = m(2) * rSges(2,1);
t1298 = rSges(3,1) * t1185 - rSges(3,2) * t1183;
t1475 = t1298 * m(3);
t1407 = ((t1205 + t1475) * t1186 - t1155 * t1184) * t1180;
t1462 = t1157 / 0.4e1;
t1463 = -t1156 / 0.4e1;
t1314 = t1180 * t1362;
t1346 = t1015 * t1447;
t1381 = t1180 * t1186;
t1441 = pkin(5) * t1031;
t1451 = pkin(2) * t1169;
t1455 = pkin(2) * t1015;
t1342 = t1183 * t1441;
t977 = t1342 - t1455;
t941 = (((t1015 * t1182 + t1031 * t1381) * t1451 - (t1346 - t1441) * t1314 + t1182 * t977) * t1031 - (-t1015 * t1381 + (-t1169 * t1182 + t1183 * t1314 + t1182) * t1031) * t1455) * t1403;
t1364 = t1182 * t1229;
t1382 = t1180 * t1183;
t1227 = pkin(5) ^ 2;
t1228 = pkin(2) ^ 2;
t1432 = t1031 * (-pkin(5) * t1346 + (t1169 * t1228 + t1227) * t1031);
t949 = (t1364 * t1432 + (-t1015 * t1133 * t1382 + t1182 * (t1015 * t1451 - t1342)) * t1015) * t1403;
t953 = (t977 * t1455 - t1432) * t1468;
t1510 = -t1057 * t941 - t1121 * t949 - t953 * t1407 - 0.4e1 * t1015 * ((t1183 * t1463 + t1185 * t1462) * t1015 + (t1318 / 0.2e1 + (t1169 - 0.1e1 / 0.2e1) * t1158) * t1031);
t1030 = t1031 ^ 2;
t1152 = t1355 * m(3) + Icges(3,3);
t1141 = rSges(3,1) * t1183 + rSges(3,2) * t1185;
t1089 = -t1180 * t1184 * t1141 + t1298 * t1182;
t1461 = m(3) * t1089;
t1466 = 0.2e1 * t1158;
t1509 = -t1121 * t941 - t1152 * t949 - t953 * t1461 + t1030 * (t1169 * t1466 - t1158 + t1318);
t1192 = sin(qJ(2,3));
t1198 = cos(qJ(2,3));
t1367 = t1182 * t1198;
t1372 = t1182 * t1192;
t1197 = cos(qJ(3,3));
t1444 = pkin(2) * t1197;
t1071 = (-t1179 * t1192 + t1181 * t1367) * t1444 + pkin(5) * (t1179 * t1198 + t1181 * t1372);
t1074 = (t1179 * t1367 + t1181 * t1192) * t1444 + (t1179 * t1372 - t1181 * t1198) * pkin(5);
t1218 = koppelP(3,2);
t1222 = koppelP(3,1);
t1126 = t1167 * t1222 + t1168 * t1218;
t1130 = -t1167 * t1218 + t1168 * t1222;
t1188 = legFrame(3,2);
t1160 = sin(t1188);
t1164 = cos(t1188);
t1281 = (t1130 * t1206 + t1208) * t1160 - (-t1126 * t1206 + t1209) * t1164;
t1360 = t1192 * t1197;
t1135 = pkin(2) * t1360 - pkin(5) * t1198;
t1191 = sin(qJ(3,3));
t1445 = pkin(2) * t1191;
t1110 = t1135 * t1180 + t1182 * t1445;
t1467 = 0.1e1 / t1110;
t1402 = t1467 / t1197;
t1325 = t1229 * t1402;
t1023 = (t1281 * t1071 + t1074 * t1207) * t1325;
t1245 = t1180 * t1197 + t1191 * t1372;
t1361 = t1191 * t1198;
t1083 = t1245 * t1179 - t1181 * t1361;
t1086 = t1179 * t1361 + t1245 * t1181;
t1038 = (t1083 * t1207 + t1281 * t1086) * t1402;
t1214 = 0.2e1 * qJ(3,3);
t1058 = cos(t1214) * t1464 - t1158 * sin(t1214) + t1230;
t1122 = -t1156 * t1197 - t1157 * t1191;
t1172 = t1197 ^ 2;
t1317 = t1147 * t1191 * t1197;
t1297 = rSges(3,1) * t1197 - rSges(3,2) * t1191;
t1476 = t1297 * m(3);
t1406 = ((t1205 + t1476) * t1198 - t1155 * t1192) * t1180;
t1313 = t1180 * t1360;
t1345 = t1023 * t1445;
t1377 = t1180 * t1198;
t1440 = pkin(5) * t1038;
t1450 = pkin(2) * t1172;
t1454 = pkin(2) * t1023;
t1341 = t1191 * t1440;
t979 = t1341 - t1454;
t942 = (((t1023 * t1182 + t1038 * t1377) * t1450 - (t1345 - t1440) * t1313 + t1182 * t979) * t1038 - (-t1023 * t1377 + (-t1172 * t1182 + t1191 * t1313 + t1182) * t1038) * t1454) * t1402;
t1380 = t1180 * t1191;
t1431 = t1038 * (-pkin(5) * t1345 + (t1172 * t1228 + t1227) * t1038);
t950 = (t1364 * t1431 + (-t1023 * t1135 * t1380 + t1182 * (t1023 * t1450 - t1341)) * t1023) * t1402;
t954 = (t979 * t1454 - t1431) * t1467;
t1508 = -t1058 * t942 - t1122 * t950 - t954 * t1406 - 0.4e1 * t1023 * ((t1191 * t1463 + t1197 * t1462) * t1023 + (t1317 / 0.2e1 + (t1172 - 0.1e1 / 0.2e1) * t1158) * t1038);
t1194 = sin(qJ(2,2));
t1200 = cos(qJ(2,2));
t1366 = t1182 * t1200;
t1370 = t1182 * t1194;
t1199 = cos(qJ(3,2));
t1443 = pkin(2) * t1199;
t1072 = (-t1179 * t1194 + t1181 * t1366) * t1443 + pkin(5) * (t1179 * t1200 + t1181 * t1370);
t1075 = (t1179 * t1366 + t1181 * t1194) * t1443 + (t1179 * t1370 - t1181 * t1200) * pkin(5);
t1219 = koppelP(2,2);
t1223 = koppelP(2,1);
t1127 = t1167 * t1223 + t1168 * t1219;
t1131 = -t1167 * t1219 + t1168 * t1223;
t1189 = legFrame(2,2);
t1161 = sin(t1189);
t1165 = cos(t1189);
t1280 = (t1131 * t1206 + t1208) * t1161 - (-t1127 * t1206 + t1209) * t1165;
t1175 = 0.1e1 / t1199;
t1358 = t1194 * t1199;
t1312 = t1180 * t1358;
t1193 = sin(qJ(3,2));
t1371 = t1182 * t1193;
t1376 = t1180 * t1200;
t1470 = 0.1e1 / (-pkin(5) * t1376 + (t1312 + t1371) * pkin(2));
t1401 = t1470 * t1175;
t1322 = t1229 * t1401;
t1024 = (t1280 * t1072 + t1075 * t1207) * t1322;
t1244 = t1180 * t1199 + t1193 * t1370;
t1359 = t1193 * t1200;
t1084 = t1244 * t1179 - t1181 * t1359;
t1087 = t1179 * t1359 + t1244 * t1181;
t1039 = (t1084 * t1207 + t1280 * t1087) * t1401;
t1215 = 0.2e1 * qJ(3,2);
t1059 = cos(t1215) * t1464 - t1158 * sin(t1215) + t1230;
t1123 = -t1156 * t1199 - t1157 * t1193;
t1174 = t1199 ^ 2;
t1316 = t1147 * t1193 * t1199;
t1296 = rSges(3,1) * t1199 - rSges(3,2) * t1193;
t1477 = t1296 * m(3);
t1405 = ((t1205 + t1477) * t1200 - t1155 * t1194) * t1180;
t1453 = pkin(2) * t1024;
t1344 = t1193 * t1453;
t1425 = t1039 * t1470;
t1439 = pkin(5) * t1039;
t1449 = pkin(2) * t1174;
t1340 = t1193 * t1439;
t980 = t1340 - t1453;
t943 = (((t1024 * t1182 + t1039 * t1376) * t1449 - (t1344 - t1439) * t1312 + t1182 * t980) * t1425 + (t1024 * t1376 + (t1174 * t1182 - t1193 * t1312 - t1182) * t1039) * t1470 * t1453) * t1175;
t1136 = pkin(2) * t1358 - pkin(5) * t1200;
t1111 = pkin(2) * t1371 + t1136 * t1180;
t1108 = 0.1e1 / t1111;
t1379 = t1180 * t1193;
t963 = -pkin(5) * t1344 + (t1174 * t1228 + t1227) * t1039;
t951 = (t963 * t1364 * t1425 + (-t1024 * t1136 * t1379 + t1182 * (t1024 * t1449 - t1340)) * t1108 * t1024) * t1175;
t955 = (t1039 * t963 - t980 * t1453) * t1470;
t1507 = -t1059 * t943 - t1123 * t951 + t955 * t1405 - 0.4e1 * t1024 * ((t1193 * t1463 + t1199 * t1462) * t1024 + (t1316 / 0.2e1 + (t1174 - 0.1e1 / 0.2e1) * t1158) * t1039);
t1196 = sin(qJ(2,1));
t1202 = cos(qJ(2,1));
t1365 = t1182 * t1202;
t1368 = t1182 * t1196;
t1201 = cos(qJ(3,1));
t1442 = pkin(2) * t1201;
t1073 = (-t1179 * t1196 + t1181 * t1365) * t1442 + pkin(5) * (t1179 * t1202 + t1181 * t1368);
t1076 = (t1179 * t1365 + t1181 * t1196) * t1442 + (t1179 * t1368 - t1181 * t1202) * pkin(5);
t1220 = koppelP(1,2);
t1224 = koppelP(1,1);
t1128 = t1167 * t1224 + t1168 * t1220;
t1132 = -t1167 * t1220 + t1168 * t1224;
t1190 = legFrame(1,2);
t1162 = sin(t1190);
t1166 = cos(t1190);
t1279 = (t1132 * t1206 + t1208) * t1162 - (-t1128 * t1206 + t1209) * t1166;
t1177 = 0.1e1 / t1201;
t1356 = t1196 * t1201;
t1311 = t1180 * t1356;
t1195 = sin(qJ(3,1));
t1369 = t1182 * t1195;
t1375 = t1180 * t1202;
t1469 = 0.1e1 / (-pkin(5) * t1375 + (t1311 + t1369) * pkin(2));
t1400 = t1469 * t1177;
t1319 = t1229 * t1400;
t1025 = (t1279 * t1073 + t1076 * t1207) * t1319;
t1243 = t1180 * t1201 + t1195 * t1368;
t1357 = t1195 * t1202;
t1085 = t1243 * t1179 - t1181 * t1357;
t1088 = t1179 * t1357 + t1243 * t1181;
t1040 = (t1085 * t1207 + t1279 * t1088) * t1400;
t1216 = 0.2e1 * qJ(3,1);
t1060 = cos(t1216) * t1464 - t1158 * sin(t1216) + t1230;
t1124 = -t1156 * t1201 - t1157 * t1195;
t1176 = t1201 ^ 2;
t1315 = t1147 * t1195 * t1201;
t1295 = rSges(3,1) * t1201 - rSges(3,2) * t1195;
t1478 = t1295 * m(3);
t1404 = ((t1205 + t1478) * t1202 - t1155 * t1196) * t1180;
t1452 = pkin(2) * t1025;
t1343 = t1195 * t1452;
t1424 = t1040 * t1469;
t1438 = pkin(5) * t1040;
t1448 = pkin(2) * t1176;
t1339 = t1195 * t1438;
t981 = t1339 - t1452;
t944 = (((t1025 * t1182 + t1040 * t1375) * t1448 - (t1343 - t1438) * t1311 + t1182 * t981) * t1424 + (t1025 * t1375 + (t1176 * t1182 - t1195 * t1311 - t1182) * t1040) * t1469 * t1452) * t1177;
t1137 = pkin(2) * t1356 - pkin(5) * t1202;
t1112 = pkin(2) * t1369 + t1137 * t1180;
t1109 = 0.1e1 / t1112;
t1378 = t1180 * t1195;
t964 = -pkin(5) * t1343 + (t1176 * t1228 + t1227) * t1040;
t952 = (t964 * t1364 * t1424 + (-t1025 * t1137 * t1378 + t1182 * (t1025 * t1448 - t1339)) * t1109 * t1025) * t1177;
t956 = (t1040 * t964 - t981 * t1452) * t1469;
t1506 = -t1060 * t944 - t1124 * t952 + t956 * t1404 - 0.4e1 * t1025 * ((t1195 * t1463 + t1201 * t1462) * t1025 + (t1315 / 0.2e1 + (t1176 - 0.1e1 / 0.2e1) * t1158) * t1040);
t1035 = t1038 ^ 2;
t1142 = rSges(3,1) * t1191 + rSges(3,2) * t1197;
t1090 = -t1180 * t1192 * t1142 + t1297 * t1182;
t1460 = m(3) * t1090;
t1505 = -t1122 * t942 - t1152 * t950 - t954 * t1460 + t1035 * (t1172 * t1466 - t1158 + t1317);
t1036 = t1039 ^ 2;
t1143 = rSges(3,1) * t1193 + rSges(3,2) * t1199;
t1091 = -t1180 * t1194 * t1143 + t1296 * t1182;
t1459 = m(3) * t1091;
t1504 = -t1123 * t943 - t1152 * t951 + t955 * t1459 + t1036 * (t1174 * t1466 - t1158 + t1316);
t1037 = t1040 ^ 2;
t1144 = rSges(3,1) * t1195 + rSges(3,2) * t1201;
t1092 = -t1180 * t1196 * t1144 + t1295 * t1182;
t1458 = m(3) * t1092;
t1503 = -t1124 * t944 - t1152 * t952 + t956 * t1458 + t1037 * (t1176 * t1466 - t1158 + t1315);
t1498 = t1503 * t1319;
t1497 = t1504 * t1322;
t1496 = t1505 * t1325;
t1495 = t1509 * t1328;
t1494 = t1506 * t1400;
t1493 = t1507 * t1401;
t1492 = t1508 * t1402;
t1491 = t1510 * t1403;
t1014 = t1015 ^ 2;
t1171 = m(1) + m(2) + m(3);
t1430 = 2 * m(3);
t1457 = m(3) * t1182;
t1490 = -t1171 * t953 - t941 * t1407 - t949 * t1461 + ((-t1030 * t1205 - (t1030 + t1014) * t1475) * t1184 - t1186 * t1031 * (t1141 * t1015 * t1430 + t1031 * t1155)) * t1180 - t1014 * t1141 * t1457;
t1020 = t1023 ^ 2;
t1489 = -t1171 * t954 - t942 * t1406 - t950 * t1460 + ((-t1035 * t1205 - (t1035 + t1020) * t1476) * t1192 - t1198 * t1038 * (t1142 * t1023 * t1430 + t1038 * t1155)) * t1180 - t1020 * t1142 * t1457;
t1021 = t1024 ^ 2;
t1488 = t1171 * t955 - t943 * t1405 - t951 * t1459 + ((-t1036 * t1205 - (t1036 + t1021) * t1477) * t1194 - t1200 * t1039 * (t1143 * t1024 * t1430 + t1039 * t1155)) * t1180 - t1021 * t1143 * t1457;
t1022 = t1025 ^ 2;
t1487 = t1171 * t956 - t944 * t1404 - t952 * t1458 + ((-t1037 * t1205 - (t1037 + t1022) * t1478) * t1196 - t1202 * t1040 * (t1144 * t1025 * t1430 + t1040 * t1155)) * t1180 - t1022 * t1144 * t1457;
t1486 = t1108 * t1488;
t1485 = t1109 * t1487;
t1484 = t1467 * t1489;
t1483 = t1468 * t1490;
t1482 = -t1498 * t1073 - t1494 * t1088;
t1481 = -t1497 * t1072 - t1493 * t1087;
t1480 = -t1496 * t1071 - t1492 * t1086;
t1479 = -t1495 * t1069 - t1491 * t1079;
t1474 = (-t1128 * t1162 + t1132 * t1166) * t1400;
t1473 = (-t1127 * t1161 + t1131 * t1165) * t1401;
t1472 = (-t1126 * t1160 + t1130 * t1164) * t1402;
t1471 = (-t1125 * t1159 + t1129 * t1163) * t1403;
t1178 = t1206 ^ 2;
t1456 = m(3) * t1229;
t1437 = m(4) * t1178;
t1411 = t1069 * t1229;
t1410 = t1071 * t1229;
t1409 = t1072 * t1229;
t1408 = t1073 * t1229;
t1399 = t1468 * t1171;
t1398 = t1467 * t1171;
t1397 = t1108 * t1171;
t1396 = t1109 * t1171;
t1395 = t1121 * t1229;
t1394 = t1122 * t1229;
t1393 = t1123 * t1229;
t1392 = t1124 * t1229;
t1383 = t1152 * t1229;
t1354 = t1468 * t1461;
t1353 = t1089 * t1456;
t1352 = t1467 * t1460;
t1351 = t1090 * t1456;
t1350 = t1108 * t1459;
t1349 = t1091 * t1456;
t1348 = t1109 * t1458;
t1347 = t1092 * t1456;
t1334 = t1468 * t1407;
t1333 = t1467 * t1406;
t1332 = t1108 * t1405;
t1331 = t1109 * t1404;
t1330 = t1159 * t1403;
t1329 = t1163 * t1403;
t1327 = t1160 * t1402;
t1326 = t1164 * t1402;
t1324 = t1161 * t1401;
t1323 = t1165 * t1401;
t1321 = t1162 * t1400;
t1320 = t1166 * t1400;
t1278 = pkin(2) * t1382 - t1133 * t1182;
t1277 = pkin(2) * t1380 - t1135 * t1182;
t1276 = pkin(2) * t1379 - t1136 * t1182;
t1275 = pkin(2) * t1378 - t1137 * t1182;
t1254 = t1057 * t1079 + t1069 * t1395;
t1253 = t1058 * t1086 + t1071 * t1394;
t1252 = t1059 * t1087 + t1072 * t1393;
t1251 = t1060 * t1088 + t1073 * t1392;
t1250 = t1069 * t1383 + t1079 * t1121;
t1249 = t1071 * t1383 + t1086 * t1122;
t1248 = t1072 * t1383 + t1087 * t1123;
t1247 = t1073 * t1383 + t1088 * t1124;
t1134 = pkin(5) * t1184 + t1186 * t1446;
t1061 = t1134 * t1181 + t1278 * t1179;
t1049 = t1061 * t1163 + t1106 * t1159;
t1050 = -t1061 * t1159 + t1106 * t1163;
t1242 = t1468 * (-t1049 * t1129 - t1050 * t1125);
t1138 = pkin(5) * t1192 + t1198 * t1444;
t1063 = t1138 * t1181 + t1277 * t1179;
t1051 = t1063 * t1164 + t1110 * t1160;
t1052 = -t1063 * t1160 + t1110 * t1164;
t1241 = t1467 * (-t1051 * t1130 - t1052 * t1126);
t1139 = pkin(5) * t1194 + t1200 * t1443;
t1064 = t1139 * t1181 + t1276 * t1179;
t1053 = t1064 * t1165 + t1111 * t1161;
t1054 = -t1064 * t1161 + t1111 * t1165;
t1240 = t1108 * (-t1053 * t1131 - t1054 * t1127);
t1140 = pkin(5) * t1196 + t1202 * t1442;
t1065 = t1140 * t1181 + t1275 * t1179;
t1055 = t1065 * t1166 + t1112 * t1162;
t1056 = -t1065 * t1162 + t1112 * t1166;
t1239 = t1109 * (-t1055 * t1132 - t1056 * t1128);
t1238 = (t1125 * t1163 + t1129 * t1159) * t1403;
t1237 = (t1126 * t1164 + t1130 * t1160) * t1402;
t1236 = (t1127 * t1165 + t1131 * t1161) * t1401;
t1235 = (t1128 * t1166 + t1132 * t1162) * t1400;
t1234 = t1069 * t1353 + t1079 * t1407;
t1233 = t1071 * t1351 + t1086 * t1406;
t1232 = t1072 * t1349 + t1087 * t1405;
t1231 = t1073 * t1347 + t1088 * t1404;
t1212 = rSges(4,1);
t1211 = rSges(4,2);
t1068 = -t1140 * t1179 + t1275 * t1181;
t1067 = -t1139 * t1179 + t1276 * t1181;
t1066 = -t1138 * t1179 + t1277 * t1181;
t1062 = -t1134 * t1179 + t1278 * t1181;
t1048 = t1088 * t1235;
t1047 = t1087 * t1236;
t1046 = t1086 * t1237;
t1045 = t1079 * t1238;
t1044 = t1235 * t1408;
t1043 = t1236 * t1409;
t1042 = t1237 * t1410;
t1041 = t1238 * t1411;
t1034 = (-t1055 * t1128 + t1056 * t1132) * t1109;
t1033 = (-t1053 * t1127 + t1054 * t1131) * t1108;
t1032 = (-t1051 * t1126 + t1052 * t1130) * t1467;
t1029 = (-t1049 * t1125 + t1050 * t1129) * t1468;
t1 = [(t1167 * t1211 - t1168 * t1212) * t1437 + t1049 * t1483 + t1051 * t1484 + t1053 * t1486 + t1055 * t1485 + ((t1051 * t1398 - t1233 * t1326) * t1241 + ((t1051 * t1352 - t1249 * t1326) * t1410 + t1086 * (t1051 * t1333 - t1253 * t1326)) * t1472 + (t1053 * t1397 - t1232 * t1323) * t1240 + ((t1053 * t1350 - t1248 * t1323) * t1409 + t1087 * (t1053 * t1332 - t1252 * t1323)) * t1473 + (t1055 * t1396 - t1231 * t1320) * t1239 + ((t1055 * t1348 - t1247 * t1320) * t1408 + t1088 * (t1055 * t1331 - t1251 * t1320)) * t1474 + (t1049 * t1399 - t1234 * t1329) * t1242 + (t1079 * (t1049 * t1334 - t1254 * t1329) + (t1049 * t1354 - t1250 * t1329) * t1411) * t1471) * t1178 + t1482 * t1166 + t1481 * t1165 + t1480 * t1164 + t1479 * t1163; -(t1167 * t1212 + t1168 * t1211) * t1437 + t1050 * t1483 + t1052 * t1484 + t1054 * t1486 + t1056 * t1485 + ((t1054 * t1397 + t1232 * t1324) * t1240 + ((t1054 * t1350 + t1248 * t1324) * t1409 + t1087 * (t1054 * t1332 + t1252 * t1324)) * t1473 + (t1056 * t1396 + t1231 * t1321) * t1239 + ((t1056 * t1348 + t1247 * t1321) * t1408 + t1088 * (t1056 * t1331 + t1251 * t1321)) * t1474 + (t1050 * t1399 + t1234 * t1330) * t1242 + (t1079 * (t1050 * t1334 + t1254 * t1330) + (t1050 * t1354 + t1250 * t1330) * t1411) * t1471 + (t1052 * t1398 + t1233 * t1327) * t1241 + ((t1052 * t1352 + t1249 * t1327) * t1410 + t1086 * (t1052 * t1333 + t1253 * t1327)) * t1472) * t1178 - t1482 * t1162 - t1481 * t1161 - t1480 * t1160 - t1479 * t1159; t1062 * t1483 + t1066 * t1484 + t1067 * t1486 + t1068 * t1485 + t1494 * t1085 + t1493 * t1084 + t1492 * t1083 + t1491 * t1078 + t1498 * t1076 + t1497 * t1075 + t1496 * t1074 + t1495 * t1070 + ((t1067 * t1397 + (t1075 * t1349 + t1084 * t1405) * t1401) * t1240 + ((t1067 * t1332 + (t1059 * t1084 + t1075 * t1393) * t1401) * t1087 + (t1067 * t1350 + (t1075 * t1383 + t1084 * t1123) * t1401) * t1409) * t1473 + (t1068 * t1396 + (t1076 * t1347 + t1085 * t1404) * t1400) * t1239 + ((t1068 * t1331 + (t1060 * t1085 + t1076 * t1392) * t1400) * t1088 + (t1068 * t1348 + (t1076 * t1383 + t1085 * t1124) * t1400) * t1408) * t1474 + (t1062 * t1399 + (t1070 * t1353 + t1078 * t1407) * t1403) * t1242 + ((t1062 * t1334 + (t1057 * t1078 + t1070 * t1395) * t1403) * t1079 + (t1062 * t1354 + (t1070 * t1383 + t1078 * t1121) * t1403) * t1411) * t1471 + (t1066 * t1398 + (t1074 * t1351 + t1083 * t1406) * t1402) * t1241 + ((t1066 * t1333 + (t1058 * t1083 + t1074 * t1394) * t1402) * t1086 + (t1066 * t1352 + (t1074 * t1383 + t1083 * t1122) * t1402) * t1410) * t1472) * t1178; t1506 * t1048 + t1507 * t1047 + t1508 * t1046 + t1510 * t1045 + t1503 * t1044 + t1504 * t1043 + t1505 * t1042 + t1509 * t1041 + t1487 * t1034 + t1488 * t1033 + t1489 * t1032 + t1490 * t1029 + ((t1033 * t1171 + t1043 * t1459 + t1047 * t1405) * t1240 + (t1087 * (t1033 * t1405 + t1043 * t1123 + t1047 * t1059) + (t1033 * t1459 + t1043 * t1152 + t1047 * t1123) * t1409) * t1473 + (t1034 * t1171 + t1044 * t1458 + t1048 * t1404) * t1239 + (t1088 * (t1034 * t1404 + t1044 * t1124 + t1048 * t1060) + (t1034 * t1458 + t1044 * t1152 + t1048 * t1124) * t1408) * t1474 + (t1029 * t1171 + t1041 * t1461 + t1045 * t1407) * t1242 + (t1079 * (t1029 * t1407 + t1041 * t1121 + t1045 * t1057) + (t1029 * t1461 + t1041 * t1152 + t1045 * t1121) * t1411) * t1471 + (t1032 * t1171 + t1042 * t1460 + t1046 * t1406) * t1241 + (t1086 * (t1032 * t1406 + t1042 * t1122 + t1046 * t1058) + (t1032 * t1460 + t1042 * t1152 + t1046 * t1122) * t1410) * t1472) * t1178;];
taucX  = t1;
