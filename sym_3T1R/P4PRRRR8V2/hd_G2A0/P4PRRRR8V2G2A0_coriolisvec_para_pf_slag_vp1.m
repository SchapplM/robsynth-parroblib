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
% Datum: 2020-08-07 11:22
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:16:14
% EndTime: 2020-08-07 11:16:33
% DurationCPUTime: 19.94s
% Computational Cost: add. (75776->631), mult. (164538->1218), div. (9024->9), fcn. (155268->30), ass. (0->464)
t1188 = sin(qJ(2,4));
t1190 = cos(qJ(2,4));
t1214 = pkin(7) + pkin(6);
t1118 = pkin(2) * t1188 - t1190 * t1214;
t1184 = sin(pkin(4));
t1186 = cos(pkin(4));
t1187 = sin(qJ(3,4));
t1333 = t1186 * t1187;
t1081 = pkin(3) * t1333 + t1118 * t1184;
t1189 = cos(qJ(3,4));
t1353 = t1184 * t1188;
t1177 = t1189 ^ 2;
t1437 = pkin(3) * t1177;
t1057 = 0.1e1 / (pkin(2) * t1333 + t1081 * t1189 + t1353 * t1437);
t1183 = sin(pkin(8));
t1185 = cos(pkin(8));
t1332 = t1186 * t1188;
t1101 = t1183 * t1332 - t1185 * t1190;
t1352 = t1184 * t1189;
t1070 = t1101 * t1187 + t1183 * t1352;
t1102 = t1183 * t1190 + t1185 * t1332;
t1337 = t1185 * t1189;
t1071 = -t1102 * t1187 - t1184 * t1337;
t1211 = xDP(3);
t1216 = xP(4);
t1175 = sin(t1216);
t1176 = cos(t1216);
t1219 = koppelP(4,2);
t1223 = koppelP(4,1);
t1109 = t1175 * t1223 + t1176 * t1219;
t1113 = -t1175 * t1219 + t1176 * t1223;
t1191 = legFrame(4,2);
t1162 = sin(t1191);
t1166 = cos(t1191);
t1210 = xDP(4);
t1212 = xDP(2);
t1213 = xDP(1);
t1255 = (t1113 * t1210 + t1212) * t1162 - (-t1109 * t1210 + t1213) * t1166;
t1014 = (t1070 * t1255 + t1071 * t1211) * t1057;
t1013 = t1014 ^ 2;
t1207 = pkin(6) + rSges(3,3);
t1447 = m(3) * t1207;
t1144 = rSges(3,2) * t1447 - Icges(3,6);
t1145 = rSges(3,1) * t1447 - Icges(3,5);
t1085 = -t1144 * t1189 - t1187 * t1145;
t1227 = rSges(3,2) ^ 2;
t1228 = rSges(3,1) ^ 2;
t1117 = (t1228 / 0.2e1 - t1227 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t1139 = (t1227 + t1228) * m(3) + Icges(3,3);
t1215 = pkin(2) * m(3);
t1310 = t1215 / 0.2e1;
t1160 = rSges(3,2) * t1310;
t1457 = m(3) * rSges(3,1);
t1311 = rSges(3,2) * t1457;
t1312 = -t1311 / 0.2e1 + Icges(3,4) / 0.2e1;
t1131 = rSges(3,1) * t1187 + rSges(3,2) * t1189;
t1263 = rSges(3,1) * t1189 - rSges(3,2) * t1187;
t1069 = -t1131 * t1353 + t1186 * t1263;
t1452 = m(3) * t1069;
t1157 = -Icges(3,4) + t1311;
t1277 = rSges(3,1) * t1310;
t1467 = t1157 * t1177 + t1187 * t1277;
t1119 = pkin(2) * t1190 + t1188 * t1214;
t1331 = t1186 * t1190;
t1338 = t1185 * t1186;
t1433 = pkin(3) * t1189;
t1045 = (t1183 * t1188 - t1185 * t1331) * t1433 - t1119 * t1338 + t1118 * t1183;
t1356 = t1183 * t1186;
t1046 = (t1183 * t1331 + t1185 * t1188) * t1433 + t1119 * t1356 + t1118 * t1185;
t1230 = 0.1e1 / pkin(3);
t1384 = t1057 * t1230;
t1006 = (t1045 * t1211 + t1046 * t1255) * t1384;
t1441 = pkin(3) * t1006;
t1305 = t1187 * t1441;
t1320 = t1188 * t1189;
t1351 = t1184 * t1190;
t1354 = t1184 * t1187;
t1411 = t1014 * t1214;
t1412 = t1014 * t1190;
t1417 = t1006 * t1057;
t1289 = t1187 * t1411;
t957 = t1289 - t1441;
t921 = ((t1006 * t1186 + t1014 * t1351) * t1437 + ((-t1305 + t1411) * t1188 + pkin(2) * t1412) * t1352 + t957 * t1186) * t1057 * t1014 + pkin(3) * (t1006 * t1351 + (t1177 * t1186 - t1320 * t1354 - t1186) * t1014) * t1417;
t1097 = pkin(3) * t1320 + t1118;
t1321 = t1186 * t1230;
t1231 = pkin(2) ^ 2;
t1158 = t1214 ^ 2 + t1231;
t1229 = pkin(3) ^ 2;
t1459 = 0.2e1 * pkin(2);
t1428 = pkin(3) * t1459;
t1427 = t1014 * (-t1214 * t1305 + (t1177 * t1229 + t1189 * t1428 + t1158) * t1014);
t925 = t1057 * t1321 * t1427 + (-t1186 * t1289 + (-t1097 * t1354 + t1186 * (pkin(2) * t1189 + t1437)) * t1006) / (t1097 * t1352 + (pkin(2) + t1433) * t1333) * t1006;
t933 = (-t1189 * t1427 - (pkin(2) * t1006 - t1189 * t957) * t1441) * t1057;
t1512 = -t1085 * t921 - t1139 * t925 - t1452 * t933 + 0.2e1 * t1013 * ((t1117 * t1187 + t1160) * t1189 + t1312 + t1467);
t1196 = sin(qJ(2,3));
t1202 = cos(qJ(2,3));
t1120 = pkin(2) * t1196 - t1202 * t1214;
t1195 = sin(qJ(3,3));
t1330 = t1186 * t1195;
t1082 = pkin(3) * t1330 + t1120 * t1184;
t1201 = cos(qJ(3,3));
t1349 = t1184 * t1196;
t1179 = t1201 ^ 2;
t1436 = pkin(3) * t1179;
t1058 = 0.1e1 / (pkin(2) * t1330 + t1082 * t1201 + t1349 * t1436);
t1329 = t1186 * t1196;
t1103 = t1183 * t1329 - t1185 * t1202;
t1344 = t1184 * t1201;
t1072 = t1103 * t1195 + t1183 * t1344;
t1106 = t1183 * t1202 + t1185 * t1329;
t1336 = t1185 * t1201;
t1075 = -t1106 * t1195 - t1184 * t1336;
t1220 = koppelP(3,2);
t1224 = koppelP(3,1);
t1110 = t1175 * t1224 + t1176 * t1220;
t1114 = -t1175 * t1220 + t1176 * t1224;
t1192 = legFrame(3,2);
t1163 = sin(t1192);
t1167 = cos(t1192);
t1254 = (t1114 * t1210 + t1212) * t1163 - (-t1110 * t1210 + t1213) * t1167;
t1018 = (t1072 * t1254 + t1075 * t1211) * t1058;
t1015 = t1018 ^ 2;
t1086 = -t1144 * t1201 - t1195 * t1145;
t1135 = rSges(3,1) * t1195 + rSges(3,2) * t1201;
t1262 = rSges(3,1) * t1201 - rSges(3,2) * t1195;
t1078 = -t1135 * t1349 + t1186 * t1262;
t1451 = m(3) * t1078;
t1466 = t1157 * t1179 + t1195 * t1277;
t1123 = pkin(2) * t1202 + t1196 * t1214;
t1324 = t1186 * t1202;
t1432 = pkin(3) * t1201;
t1047 = (t1183 * t1196 - t1185 * t1324) * t1432 - t1123 * t1338 + t1120 * t1183;
t1050 = (t1183 * t1324 + t1185 * t1196) * t1432 + t1123 * t1356 + t1120 * t1185;
t1382 = t1058 * t1230;
t1010 = (t1047 * t1211 + t1050 * t1254) * t1382;
t1440 = pkin(3) * t1010;
t1304 = t1195 * t1440;
t1319 = t1196 * t1201;
t1343 = t1184 * t1202;
t1350 = t1184 * t1195;
t1406 = t1018 * t1214;
t1407 = t1018 * t1202;
t1416 = t1010 * t1058;
t1288 = t1195 * t1406;
t959 = t1288 - t1440;
t922 = ((t1010 * t1186 + t1018 * t1343) * t1436 + ((-t1304 + t1406) * t1196 + pkin(2) * t1407) * t1344 + t959 * t1186) * t1058 * t1018 + pkin(3) * (t1010 * t1343 + (t1179 * t1186 - t1319 * t1350 - t1186) * t1018) * t1416;
t1098 = pkin(3) * t1319 + t1120;
t1426 = t1018 * (-t1214 * t1304 + (t1179 * t1229 + t1201 * t1428 + t1158) * t1018);
t926 = t1058 * t1321 * t1426 + (-t1186 * t1288 + (-t1098 * t1350 + t1186 * (pkin(2) * t1201 + t1436)) * t1010) / (t1098 * t1344 + (pkin(2) + t1432) * t1330) * t1010;
t934 = (-t1201 * t1426 - (pkin(2) * t1010 - t1201 * t959) * t1440) * t1058;
t1511 = -t1086 * t922 - t1139 * t926 - t1451 * t934 + 0.2e1 * t1015 * ((t1117 * t1195 + t1160) * t1201 + t1312 + t1466);
t1198 = sin(qJ(2,2));
t1204 = cos(qJ(2,2));
t1121 = pkin(2) * t1198 - t1204 * t1214;
t1197 = sin(qJ(3,2));
t1328 = t1186 * t1197;
t1083 = pkin(3) * t1328 + t1121 * t1184;
t1203 = cos(qJ(3,2));
t1347 = t1184 * t1198;
t1180 = t1203 ^ 2;
t1435 = pkin(3) * t1180;
t1059 = 0.1e1 / (pkin(2) * t1328 + t1083 * t1203 + t1347 * t1435);
t1327 = t1186 * t1198;
t1104 = t1183 * t1327 - t1185 * t1204;
t1342 = t1184 * t1203;
t1073 = t1104 * t1197 + t1183 * t1342;
t1107 = t1183 * t1204 + t1185 * t1327;
t1335 = t1185 * t1203;
t1076 = -t1107 * t1197 - t1184 * t1335;
t1221 = koppelP(2,2);
t1225 = koppelP(2,1);
t1111 = t1175 * t1225 + t1176 * t1221;
t1115 = -t1175 * t1221 + t1176 * t1225;
t1193 = legFrame(2,2);
t1164 = sin(t1193);
t1168 = cos(t1193);
t1253 = (t1115 * t1210 + t1212) * t1164 - (-t1111 * t1210 + t1213) * t1168;
t1019 = (t1073 * t1253 + t1076 * t1211) * t1059;
t1016 = t1019 ^ 2;
t1087 = -t1144 * t1203 - t1197 * t1145;
t1136 = rSges(3,1) * t1197 + rSges(3,2) * t1203;
t1261 = rSges(3,1) * t1203 - rSges(3,2) * t1197;
t1079 = -t1136 * t1347 + t1186 * t1261;
t1450 = m(3) * t1079;
t1465 = t1157 * t1180 + t1197 * t1277;
t1124 = pkin(2) * t1204 + t1198 * t1214;
t1323 = t1186 * t1204;
t1431 = pkin(3) * t1203;
t1048 = (t1183 * t1198 - t1185 * t1323) * t1431 - t1124 * t1338 + t1121 * t1183;
t1051 = (t1183 * t1323 + t1185 * t1198) * t1431 + t1124 * t1356 + t1121 * t1185;
t1380 = t1059 * t1230;
t1011 = (t1048 * t1211 + t1051 * t1253) * t1380;
t1439 = pkin(3) * t1011;
t1303 = t1197 * t1439;
t1318 = t1198 * t1203;
t1341 = t1184 * t1204;
t1348 = t1184 * t1197;
t1404 = t1019 * t1214;
t1405 = t1019 * t1204;
t1415 = t1011 * t1059;
t1287 = t1197 * t1404;
t960 = t1287 - t1439;
t923 = ((t1011 * t1186 + t1019 * t1341) * t1435 + ((-t1303 + t1404) * t1198 + pkin(2) * t1405) * t1342 + t960 * t1186) * t1059 * t1019 + pkin(3) * (t1011 * t1341 + (t1180 * t1186 - t1318 * t1348 - t1186) * t1019) * t1415;
t1099 = pkin(3) * t1318 + t1121;
t1425 = t1019 * (-t1214 * t1303 + (t1180 * t1229 + t1203 * t1428 + t1158) * t1019);
t927 = t1059 * t1321 * t1425 + (-t1186 * t1287 + (-t1099 * t1348 + t1186 * (pkin(2) * t1203 + t1435)) * t1011) / (t1099 * t1342 + (pkin(2) + t1431) * t1328) * t1011;
t935 = (-t1203 * t1425 - (pkin(2) * t1011 - t1203 * t960) * t1439) * t1059;
t1510 = -t1087 * t923 - t1139 * t927 - t1450 * t935 + 0.2e1 * t1016 * ((t1117 * t1197 + t1160) * t1203 + t1312 + t1465);
t1200 = sin(qJ(2,1));
t1206 = cos(qJ(2,1));
t1122 = pkin(2) * t1200 - t1206 * t1214;
t1199 = sin(qJ(3,1));
t1326 = t1186 * t1199;
t1084 = pkin(3) * t1326 + t1122 * t1184;
t1205 = cos(qJ(3,1));
t1345 = t1184 * t1200;
t1181 = t1205 ^ 2;
t1434 = pkin(3) * t1181;
t1060 = 0.1e1 / (pkin(2) * t1326 + t1084 * t1205 + t1345 * t1434);
t1325 = t1186 * t1200;
t1105 = t1183 * t1325 - t1185 * t1206;
t1340 = t1184 * t1205;
t1074 = t1105 * t1199 + t1183 * t1340;
t1108 = t1183 * t1206 + t1185 * t1325;
t1334 = t1185 * t1205;
t1077 = -t1108 * t1199 - t1184 * t1334;
t1222 = koppelP(1,2);
t1226 = koppelP(1,1);
t1112 = t1175 * t1226 + t1176 * t1222;
t1116 = -t1175 * t1222 + t1176 * t1226;
t1194 = legFrame(1,2);
t1165 = sin(t1194);
t1169 = cos(t1194);
t1252 = (t1116 * t1210 + t1212) * t1165 - (-t1112 * t1210 + t1213) * t1169;
t1020 = (t1074 * t1252 + t1077 * t1211) * t1060;
t1017 = t1020 ^ 2;
t1088 = -t1144 * t1205 - t1199 * t1145;
t1137 = rSges(3,1) * t1199 + rSges(3,2) * t1205;
t1260 = rSges(3,1) * t1205 - rSges(3,2) * t1199;
t1080 = -t1137 * t1345 + t1186 * t1260;
t1449 = m(3) * t1080;
t1464 = t1157 * t1181 + t1199 * t1277;
t1125 = pkin(2) * t1206 + t1200 * t1214;
t1322 = t1186 * t1206;
t1430 = pkin(3) * t1205;
t1049 = (t1183 * t1200 - t1185 * t1322) * t1430 - t1125 * t1338 + t1122 * t1183;
t1052 = (t1183 * t1322 + t1185 * t1200) * t1430 + t1125 * t1356 + t1122 * t1185;
t1378 = t1060 * t1230;
t1012 = (t1049 * t1211 + t1052 * t1252) * t1378;
t1438 = pkin(3) * t1012;
t1302 = t1199 * t1438;
t1317 = t1200 * t1205;
t1339 = t1184 * t1206;
t1346 = t1184 * t1199;
t1402 = t1020 * t1214;
t1403 = t1020 * t1206;
t1414 = t1012 * t1060;
t1286 = t1199 * t1402;
t961 = t1286 - t1438;
t924 = ((t1012 * t1186 + t1020 * t1339) * t1434 + ((-t1302 + t1402) * t1200 + pkin(2) * t1403) * t1340 + t961 * t1186) * t1060 * t1020 + pkin(3) * (t1012 * t1339 + (t1181 * t1186 - t1317 * t1346 - t1186) * t1020) * t1414;
t1100 = pkin(3) * t1317 + t1122;
t1424 = t1020 * (-t1214 * t1302 + (t1181 * t1229 + t1205 * t1428 + t1158) * t1020);
t928 = t1060 * t1321 * t1424 + (-t1186 * t1286 + (-t1100 * t1346 + t1186 * (pkin(2) * t1205 + t1434)) * t1012) / (t1100 * t1340 + (pkin(2) + t1430) * t1326) * t1012;
t936 = (-t1205 * t1424 - (pkin(2) * t1012 - t1205 * t961) * t1438) * t1060;
t1509 = -t1088 * t924 - t1139 * t928 - t1449 * t936 + 0.2e1 * t1017 * ((t1117 * t1199 + t1160) * t1205 + t1312 + t1464);
t1453 = -t1157 / 0.2e1;
t1454 = t1145 / 0.4e1;
t1455 = -t1144 / 0.4e1;
t1138 = (-t1227 + t1228) * m(3) + Icges(3,2) - Icges(3,1);
t1456 = t1138 / 0.2e1;
t1507 = -0.4e1 * (t1187 * t1455 + t1189 * t1454) * t1006 - 0.4e1 * ((t1187 * t1456 + t1160) * t1189 + t1453 + t1467) * t1014;
t1506 = -0.4e1 * (t1195 * t1455 + t1201 * t1454) * t1010 - 0.4e1 * ((t1195 * t1456 + t1160) * t1201 + t1453 + t1466) * t1018;
t1505 = -0.4e1 * (t1199 * t1455 + t1205 * t1454) * t1012 - 0.4e1 * ((t1199 * t1456 + t1160) * t1205 + t1453 + t1464) * t1020;
t1501 = t1510 * t1380;
t1500 = t1509 * t1378;
t1141 = t1207 ^ 2 + t1227 + t1231;
t1171 = t1457 * t1459;
t1268 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) + Icges(2,3);
t1419 = -0.2e1 * rSges(3,2) * pkin(2);
t1458 = -0.2e1 * t1157;
t1043 = t1138 * t1180 + (t1197 * t1458 + t1171) * t1203 + (t1197 * t1419 + t1141) * m(3) + t1268;
t1140 = m(2) * rSges(2,2) - t1447;
t1153 = m(2) * rSges(2,1) + t1215;
t1482 = m(3) * t1261;
t1375 = ((t1153 + t1482) * t1204 - t1198 * t1140) * t1184;
t913 = -t1043 * t923 - t1087 * t927 - t1375 * t935;
t939 = (t1197 * t1455 + t1203 * t1454) * t1011 + ((t1197 * t1456 + t1160) * t1203 + t1453 + t1465) * t1019;
t1499 = -t1059 * t913 + 0.4e1 * t939 * t1415;
t1041 = t1138 * t1177 + (t1187 * t1458 + t1171) * t1189 + (t1187 * t1419 + t1141) * m(3) + t1268;
t1480 = m(3) * t1263;
t1377 = ((t1153 + t1480) * t1190 - t1188 * t1140) * t1184;
t909 = -t1041 * t921 - t1085 * t925 - t1377 * t933;
t1498 = t1057 * t909 + t1417 * t1507;
t1042 = t1138 * t1179 + (t1195 * t1458 + t1171) * t1201 + (t1195 * t1419 + t1141) * m(3) + t1268;
t1481 = m(3) * t1262;
t1376 = ((t1153 + t1481) * t1202 - t1196 * t1140) * t1184;
t912 = -t1042 * t922 - t1086 * t926 - t1376 * t934;
t1497 = t1058 * t912 + t1416 * t1506;
t1044 = t1138 * t1181 + (t1199 * t1458 + t1171) * t1205 + (t1199 * t1419 + t1141) * m(3) + t1268;
t1483 = m(3) * t1260;
t1374 = ((t1153 + t1483) * t1206 - t1200 * t1140) * t1184;
t914 = -t1044 * t924 - t1088 * t928 - t1374 * t936;
t1496 = t1060 * t914 + t1414 * t1505;
t1495 = t1512 * t1384;
t1494 = t1511 * t1382;
t1002 = t1006 ^ 2;
t1178 = m(1) + m(2) + m(3);
t1418 = 0.2e1 * m(3);
t1448 = m(3) * t1186;
t1479 = -t1178 * t933 - t1377 * t921 - t1452 * t925 + ((-t1013 * t1153 - (t1013 + t1002) * t1480) * t1188 - (t1006 * t1131 * t1418 + t1014 * t1140) * t1412) * t1184 - t1002 * t1131 * t1448;
t1007 = t1010 ^ 2;
t1478 = -t1178 * t934 - t1376 * t922 - t1451 * t926 + ((-t1015 * t1153 - (t1015 + t1007) * t1481) * t1196 - (t1010 * t1135 * t1418 + t1018 * t1140) * t1407) * t1184 - t1007 * t1135 * t1448;
t1008 = t1011 ^ 2;
t1477 = -t1178 * t935 - t1375 * t923 - t1450 * t927 + ((-t1016 * t1153 - (t1016 + t1008) * t1482) * t1198 - (t1011 * t1136 * t1418 + t1019 * t1140) * t1405) * t1184 - t1008 * t1136 * t1448;
t1009 = t1012 ^ 2;
t1476 = -t1178 * t936 - t1374 * t924 - t1449 * t928 + ((-t1017 * t1153 - (t1017 + t1009) * t1483) * t1200 - (t1012 * t1137 * t1418 + t1020 * t1140) * t1403) * t1184 - t1009 * t1137 * t1448;
t1471 = t1057 * t1479;
t1470 = t1058 * t1478;
t1469 = t1059 * t1477;
t1468 = t1060 * t1476;
t1463 = t1500 * t1052 + t1496 * t1074;
t1462 = t1494 * t1050 + t1497 * t1072;
t1461 = -t1501 * t1051 + t1499 * t1073;
t1460 = t1495 * t1046 + t1498 * t1070;
t1182 = t1210 ^ 2;
t1446 = m(3) * t1230;
t1445 = pkin(2) * t1187;
t1444 = pkin(2) * t1195;
t1443 = pkin(2) * t1197;
t1442 = pkin(2) * t1199;
t1429 = m(4) * t1182;
t1389 = t1046 * t1230;
t1388 = t1050 * t1230;
t1387 = t1051 * t1230;
t1386 = t1052 * t1230;
t1385 = t1057 * t1182;
t1383 = t1058 * t1182;
t1381 = t1059 * t1182;
t1379 = t1060 * t1182;
t1373 = t1070 * t1162;
t1372 = t1070 * t1166;
t1371 = t1072 * t1163;
t1370 = t1072 * t1167;
t1369 = t1073 * t1164;
t1368 = t1073 * t1168;
t1367 = t1074 * t1165;
t1366 = t1074 * t1169;
t1365 = t1085 * t1230;
t1364 = t1086 * t1230;
t1363 = t1087 * t1230;
t1362 = t1088 * t1230;
t1361 = t1139 * t1230;
t1355 = t1184 * t1185;
t1309 = t1069 * t1446;
t1308 = t1078 * t1446;
t1307 = t1079 * t1446;
t1306 = t1080 * t1446;
t1285 = t1162 * t1389;
t1284 = t1166 * t1389;
t1283 = t1163 * t1388;
t1282 = t1167 * t1388;
t1281 = t1164 * t1387;
t1280 = t1168 * t1387;
t1279 = t1165 * t1386;
t1278 = t1169 * t1386;
t1251 = pkin(3) * t1354 - t1118 * t1186;
t1250 = pkin(3) * t1350 - t1120 * t1186;
t1249 = pkin(3) * t1348 - t1121 * t1186;
t1248 = pkin(3) * t1346 - t1122 * t1186;
t1247 = t1041 * t1070 + t1046 * t1365;
t1246 = t1042 * t1072 + t1050 * t1364;
t1245 = t1043 * t1073 + t1051 * t1363;
t1244 = t1044 * t1074 + t1052 * t1362;
t1243 = t1046 * t1361 + t1070 * t1085;
t1242 = t1050 * t1361 + t1072 * t1086;
t1241 = t1051 * t1361 + t1073 * t1087;
t1240 = t1052 * t1361 + t1074 * t1088;
t1239 = t1057 * (t1109 * t1166 + t1113 * t1162);
t1238 = t1058 * (t1110 * t1167 + t1114 * t1163);
t1237 = t1059 * (t1111 * t1168 + t1115 * t1164);
t1236 = t1060 * (t1112 * t1169 + t1116 * t1165);
t1235 = t1046 * t1309 + t1070 * t1377;
t1234 = t1050 * t1308 + t1072 * t1376;
t1233 = t1051 * t1307 + t1073 * t1375;
t1232 = t1052 * t1306 + t1074 * t1374;
t1218 = rSges(4,1);
t1217 = rSges(4,2);
t1064 = -t1125 * t1183 + t1185 * t1248;
t1063 = -t1124 * t1183 + t1185 * t1249;
t1062 = -t1123 * t1183 + t1185 * t1250;
t1061 = -t1119 * t1183 + t1185 * t1251;
t1040 = -t1105 * t1434 + t1125 * t1334 + (pkin(2) * t1346 + t1205 * t1248) * t1183;
t1039 = -t1104 * t1435 + t1124 * t1335 + (pkin(2) * t1348 + t1203 * t1249) * t1183;
t1038 = -t1103 * t1436 + t1123 * t1336 + (pkin(2) * t1350 + t1201 * t1250) * t1183;
t1037 = -t1101 * t1437 + t1119 * t1337 + (pkin(2) * t1354 + t1189 * t1251) * t1183;
t1036 = -(t1108 * t1165 - t1169 * t1345) * t1434 + (t1064 * t1165 + t1084 * t1169) * t1205 + (t1165 * t1355 + t1169 * t1186) * t1442;
t1035 = -(t1107 * t1164 - t1168 * t1347) * t1435 + (t1063 * t1164 + t1083 * t1168) * t1203 + (t1164 * t1355 + t1168 * t1186) * t1443;
t1034 = -(t1106 * t1163 - t1167 * t1349) * t1436 + (t1062 * t1163 + t1082 * t1167) * t1201 + (t1163 * t1355 + t1167 * t1186) * t1444;
t1033 = (t1108 * t1169 + t1165 * t1345) * t1434 + (-t1064 * t1169 + t1084 * t1165) * t1205 + (t1165 * t1186 - t1169 * t1355) * t1442;
t1032 = (t1107 * t1168 + t1164 * t1347) * t1435 + (-t1063 * t1168 + t1083 * t1164) * t1203 + (t1164 * t1186 - t1168 * t1355) * t1443;
t1031 = (t1106 * t1167 + t1163 * t1349) * t1436 + (-t1062 * t1167 + t1082 * t1163) * t1201 + (t1163 * t1186 - t1167 * t1355) * t1444;
t1030 = -(t1102 * t1162 - t1166 * t1353) * t1437 + (t1061 * t1162 + t1081 * t1166) * t1189 + (t1162 * t1355 + t1166 * t1186) * t1445;
t1029 = (t1102 * t1166 + t1162 * t1353) * t1437 + (-t1061 * t1166 + t1081 * t1162) * t1189 + (t1162 * t1186 - t1166 * t1355) * t1445;
t1028 = t1074 * t1236;
t1027 = t1073 * t1237;
t1026 = t1072 * t1238;
t1025 = t1070 * t1239;
t1024 = t1236 * t1386;
t1023 = t1237 * t1387;
t1022 = t1238 * t1388;
t1021 = t1239 * t1389;
t1005 = (t1040 * t1449 + t1049 * t1361 + t1077 * t1088) * t1060;
t1004 = (t1039 * t1450 + t1048 * t1361 + t1076 * t1087) * t1059;
t1003 = (t1038 * t1451 + t1047 * t1361 + t1075 * t1086) * t1058;
t1001 = (t1040 * t1178 + t1049 * t1306 + t1077 * t1374) * t1060;
t1000 = (t1039 * t1178 + t1048 * t1307 + t1076 * t1375) * t1059;
t999 = (t1038 * t1178 + t1047 * t1308 + t1075 * t1376) * t1058;
t998 = (t1037 * t1452 + t1045 * t1361 + t1071 * t1085) * t1057;
t997 = (t1037 * t1178 + t1045 * t1309 + t1071 * t1377) * t1057;
t996 = (-t1033 * t1112 + t1036 * t1116) * t1060;
t995 = (-t1032 * t1111 + t1035 * t1115) * t1059;
t994 = (-t1031 * t1110 + t1034 * t1114) * t1058;
t993 = (-t1029 * t1109 + t1030 * t1113) * t1057;
t992 = (t1040 * t1374 + t1044 * t1077 + t1049 * t1362) * t1060;
t991 = (t1039 * t1375 + t1043 * t1076 + t1048 * t1363) * t1059;
t990 = (t1038 * t1376 + t1042 * t1075 + t1047 * t1364) * t1058;
t989 = (t1037 * t1377 + t1041 * t1071 + t1045 * t1365) * t1057;
t988 = (t1036 * t1449 + t1165 * t1240) * t1060;
t987 = (t1035 * t1450 + t1164 * t1241) * t1059;
t986 = (t1034 * t1451 + t1163 * t1242) * t1058;
t985 = (t1033 * t1449 - t1169 * t1240) * t1060;
t984 = (t1032 * t1450 - t1168 * t1241) * t1059;
t983 = (t1031 * t1451 - t1167 * t1242) * t1058;
t982 = (t1036 * t1178 + t1165 * t1232) * t1060;
t981 = (t1035 * t1178 + t1164 * t1233) * t1059;
t980 = (t1034 * t1178 + t1163 * t1234) * t1058;
t979 = (t1033 * t1178 - t1169 * t1232) * t1060;
t978 = (t1032 * t1178 - t1168 * t1233) * t1059;
t977 = (t1031 * t1178 - t1167 * t1234) * t1058;
t976 = (t1030 * t1452 + t1162 * t1243) * t1057;
t975 = (t1029 * t1452 - t1166 * t1243) * t1057;
t974 = (t1030 * t1178 + t1162 * t1235) * t1057;
t973 = (t1029 * t1178 - t1166 * t1235) * t1057;
t972 = (t1036 * t1374 + t1165 * t1244) * t1060;
t971 = (t1035 * t1375 + t1164 * t1245) * t1059;
t970 = (t1034 * t1376 + t1163 * t1246) * t1058;
t969 = (t1033 * t1374 - t1169 * t1244) * t1060;
t968 = (t1032 * t1375 - t1168 * t1245) * t1059;
t967 = (t1031 * t1376 - t1167 * t1246) * t1058;
t966 = (t1030 * t1377 + t1162 * t1247) * t1057;
t965 = (t1029 * t1377 - t1166 * t1247) * t1057;
t956 = t1024 * t1139 + t1028 * t1088 + t1449 * t996;
t955 = t1023 * t1139 + t1027 * t1087 + t1450 * t995;
t954 = t1022 * t1139 + t1026 * t1086 + t1451 * t994;
t953 = t1024 * t1449 + t1028 * t1374 + t1178 * t996;
t952 = t1023 * t1450 + t1027 * t1375 + t1178 * t995;
t951 = t1022 * t1451 + t1026 * t1376 + t1178 * t994;
t950 = t1021 * t1139 + t1025 * t1085 + t1452 * t993;
t949 = t1021 * t1452 + t1025 * t1377 + t1178 * t993;
t948 = t1024 * t1088 + t1028 * t1044 + t1374 * t996;
t947 = t1023 * t1087 + t1027 * t1043 + t1375 * t995;
t946 = t1022 * t1086 + t1026 * t1042 + t1376 * t994;
t945 = t1021 * t1085 + t1025 * t1041 + t1377 * t993;
t1 = [(t1175 * t1217 - t1176 * t1218) * t1429 + (-(t1033 * t979 - t1278 * t985 - t1366 * t969) * t1116 - (t1036 * t979 + t1279 * t985 + t1367 * t969) * t1112) * t1379 + (-(t1032 * t978 - t1280 * t984 - t1368 * t968) * t1115 - (t1035 * t978 + t1281 * t984 + t1369 * t968) * t1111) * t1381 + (-(t1031 * t977 - t1282 * t983 - t1370 * t967) * t1114 - (t1034 * t977 + t1283 * t983 + t1371 * t967) * t1110) * t1383 + (-(t1029 * t973 - t1284 * t975 - t1372 * t965) * t1113 - (t1030 * t973 + t1285 * t975 + t1373 * t965) * t1109) * t1385 + t1029 * t1471 + t1031 * t1470 + t1032 * t1469 + t1033 * t1468 - t1463 * t1169 + t1461 * t1168 - t1462 * t1167 - t1460 * t1166; (-(t1033 * t982 - t1278 * t988 - t972 * t1366) * t1116 - (t1036 * t982 + t1279 * t988 + t972 * t1367) * t1112) * t1379 + (-(t1032 * t981 - t1280 * t987 - t1368 * t971) * t1115 - (t1035 * t981 + t1281 * t987 + t1369 * t971) * t1111) * t1381 + (-(t1031 * t980 - t1282 * t986 - t1370 * t970) * t1114 - (t1034 * t980 + t1283 * t986 + t1371 * t970) * t1110) * t1383 + (-(t1029 * t974 - t1284 * t976 - t1372 * t966) * t1113 - (t1030 * t974 + t1285 * t976 + t1373 * t966) * t1109) * t1385 - (t1175 * t1218 + t1176 * t1217) * t1429 + t1030 * t1471 + t1034 * t1470 + t1035 * t1469 + t1036 * t1468 + t1463 * t1165 - t1461 * t1164 + t1462 * t1163 + t1460 * t1162; (-(t1001 * t1033 - t1005 * t1278 - t1366 * t992) * t1116 - (t1001 * t1036 + t1005 * t1279 + t1367 * t992) * t1112) * t1379 + (-(t1000 * t1032 - t1004 * t1280 - t1368 * t991) * t1115 - (t1000 * t1035 + t1004 * t1281 + t1369 * t991) * t1111) * t1381 + (-(-t1003 * t1282 + t1031 * t999 - t1370 * t990) * t1114 - (t1003 * t1283 + t1034 * t999 + t1371 * t990) * t1110) * t1383 + (-(t1029 * t997 - t1284 * t998 - t1372 * t989) * t1113 - (t1030 * t997 + t1285 * t998 + t1373 * t989) * t1109) * t1385 + t1037 * t1471 + t1038 * t1470 + t1039 * t1469 + t1040 * t1468 + t1496 * t1077 - t1499 * t1076 + t1497 * t1075 + t1498 * t1071 + t1500 * t1049 + t1501 * t1048 + t1494 * t1047 + t1495 * t1045; (-(t1033 * t953 - t1278 * t956 - t1366 * t948) * t1116 - (t1036 * t953 + t1279 * t956 + t1367 * t948) * t1112) * t1379 + (-(t1032 * t952 - t1280 * t955 - t1368 * t947) * t1115 - (t1035 * t952 + t1281 * t955 + t1369 * t947) * t1111) * t1381 + (-(t1031 * t951 - t1282 * t954 - t1370 * t946) * t1114 - (t1034 * t951 + t1283 * t954 + t1371 * t946) * t1110) * t1383 + (-(t1029 * t949 - t1284 * t950 - t1372 * t945) * t1113 - (t1030 * t949 + t1285 * t950 + t1373 * t945) * t1109) * t1385 + t1476 * t996 + t1477 * t995 + t1478 * t994 + t1479 * t993 + (t1012 * t1505 + t914) * t1028 + (-0.4e1 * t1011 * t939 + t913) * t1027 + (t1010 * t1506 + t912) * t1026 + (t1006 * t1507 + t909) * t1025 + t1509 * t1024 + t1510 * t1023 + t1511 * t1022 + t1512 * t1021;];
taucX  = t1;
