% Calculate vector of centrifugal and coriolis load on the joints for
% P4PRRRR8V2G3A0
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
% Datum: 2020-08-07 11:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:25:11
% EndTime: 2020-08-07 11:25:28
% DurationCPUTime: 17.59s
% Computational Cost: add. (75776->643), mult. (164538->1238), div. (9024->9), fcn. (155268->30), ass. (0->472)
t1197 = sin(qJ(2,1));
t1203 = cos(qJ(2,1));
t1211 = pkin(7) + pkin(6);
t1119 = pkin(2) * t1197 - t1203 * t1211;
t1181 = sin(pkin(4));
t1183 = cos(pkin(4));
t1196 = sin(qJ(3,1));
t1323 = t1183 * t1196;
t1081 = pkin(3) * t1323 + t1119 * t1181;
t1202 = cos(qJ(3,1));
t1338 = t1181 * t1197;
t1178 = t1202 ^ 2;
t1427 = pkin(3) * t1178;
t1057 = 0.1e1 / (pkin(2) * t1323 + t1081 * t1202 + t1338 * t1427);
t1180 = sin(pkin(8));
t1182 = cos(pkin(8));
t1322 = t1183 * t1197;
t1102 = t1180 * t1322 - t1182 * t1203;
t1333 = t1181 * t1202;
t1071 = t1102 * t1196 + t1180 * t1333;
t1105 = t1180 * t1203 + t1182 * t1322;
t1074 = t1105 * t1196 + t1182 * t1333;
t1208 = xDP(3);
t1213 = xP(4);
t1172 = sin(t1213);
t1173 = cos(t1213);
t1219 = koppelP(1,2);
t1223 = koppelP(1,1);
t1109 = t1172 * t1223 + t1173 * t1219;
t1113 = -t1172 * t1219 + t1173 * t1223;
t1191 = legFrame(1,2);
t1162 = sin(t1191);
t1166 = cos(t1191);
t1207 = xDP(4);
t1209 = xDP(2);
t1210 = xDP(1);
t1249 = (t1113 * t1207 + t1209) * t1162 - (-t1109 * t1207 + t1210) * t1166;
t1017 = (t1071 * t1208 + t1074 * t1249) * t1057;
t1014 = t1017 ^ 2;
t1224 = rSges(3,2) ^ 2;
t1225 = rSges(3,1) ^ 2;
t1114 = (t1225 / 0.2e1 - t1224 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t1212 = pkin(2) * m(3);
t1307 = t1212 / 0.2e1;
t1157 = rSges(3,2) * t1307;
t1450 = m(3) * rSges(3,1);
t1308 = rSges(3,2) * t1450;
t1309 = -t1308 / 0.2e1 + Icges(3,4) / 0.2e1;
t1154 = -Icges(3,4) + t1308;
t1274 = rSges(3,1) * t1307;
t1457 = t1154 * t1178 + t1196 * t1274;
t1488 = 0.2e1 * t1014 * ((t1114 * t1196 + t1157) * t1202 + t1309 + t1457);
t1195 = sin(qJ(2,2));
t1201 = cos(qJ(2,2));
t1118 = pkin(2) * t1195 - t1201 * t1211;
t1194 = sin(qJ(3,2));
t1325 = t1183 * t1194;
t1080 = pkin(3) * t1325 + t1118 * t1181;
t1200 = cos(qJ(3,2));
t1340 = t1181 * t1195;
t1177 = t1200 ^ 2;
t1428 = pkin(3) * t1177;
t1056 = 0.1e1 / (pkin(2) * t1325 + t1080 * t1200 + t1340 * t1428);
t1324 = t1183 * t1195;
t1101 = t1180 * t1324 - t1182 * t1201;
t1335 = t1181 * t1200;
t1070 = t1101 * t1194 + t1180 * t1335;
t1104 = t1180 * t1201 + t1182 * t1324;
t1073 = t1104 * t1194 + t1182 * t1335;
t1218 = koppelP(2,2);
t1222 = koppelP(2,1);
t1108 = t1172 * t1222 + t1173 * t1218;
t1112 = -t1172 * t1218 + t1173 * t1222;
t1190 = legFrame(2,2);
t1161 = sin(t1190);
t1165 = cos(t1190);
t1250 = (t1112 * t1207 + t1209) * t1161 - (-t1108 * t1207 + t1210) * t1165;
t1016 = (t1070 * t1208 + t1073 * t1250) * t1056;
t1013 = t1016 ^ 2;
t1458 = t1154 * t1177 + t1194 * t1274;
t1487 = 0.2e1 * t1013 * ((t1114 * t1194 + t1157) * t1200 + t1309 + t1458);
t1193 = sin(qJ(2,3));
t1199 = cos(qJ(2,3));
t1117 = pkin(2) * t1193 - t1199 * t1211;
t1192 = sin(qJ(3,3));
t1327 = t1183 * t1192;
t1079 = pkin(3) * t1327 + t1117 * t1181;
t1198 = cos(qJ(3,3));
t1342 = t1181 * t1193;
t1176 = t1198 ^ 2;
t1429 = pkin(3) * t1176;
t1055 = 0.1e1 / (pkin(2) * t1327 + t1079 * t1198 + t1342 * t1429);
t1326 = t1183 * t1193;
t1100 = t1180 * t1326 - t1182 * t1199;
t1337 = t1181 * t1198;
t1069 = t1100 * t1192 + t1180 * t1337;
t1103 = t1180 * t1199 + t1182 * t1326;
t1072 = t1103 * t1192 + t1182 * t1337;
t1217 = koppelP(3,2);
t1221 = koppelP(3,1);
t1107 = t1172 * t1221 + t1173 * t1217;
t1111 = -t1172 * t1217 + t1173 * t1221;
t1189 = legFrame(3,2);
t1160 = sin(t1189);
t1164 = cos(t1189);
t1251 = (t1111 * t1207 + t1209) * t1160 - (-t1107 * t1207 + t1210) * t1164;
t1015 = (t1069 * t1208 + t1072 * t1251) * t1055;
t1012 = t1015 ^ 2;
t1459 = t1154 * t1176 + t1192 * t1274;
t1486 = 0.2e1 * t1012 * ((t1114 * t1192 + t1157) * t1198 + t1309 + t1459);
t1185 = sin(qJ(2,4));
t1187 = cos(qJ(2,4));
t1115 = pkin(2) * t1185 - t1187 * t1211;
t1184 = sin(qJ(3,4));
t1330 = t1183 * t1184;
t1078 = pkin(3) * t1330 + t1115 * t1181;
t1186 = cos(qJ(3,4));
t1346 = t1181 * t1185;
t1174 = t1186 ^ 2;
t1430 = pkin(3) * t1174;
t1054 = 0.1e1 / (pkin(2) * t1330 + t1078 * t1186 + t1346 * t1430);
t1329 = t1183 * t1185;
t1098 = t1180 * t1329 - t1182 * t1187;
t1345 = t1181 * t1186;
t1067 = t1098 * t1184 + t1180 * t1345;
t1099 = t1180 * t1187 + t1182 * t1329;
t1068 = t1099 * t1184 + t1182 * t1345;
t1216 = koppelP(4,2);
t1220 = koppelP(4,1);
t1106 = t1172 * t1220 + t1173 * t1216;
t1110 = -t1172 * t1216 + t1173 * t1220;
t1188 = legFrame(4,2);
t1159 = sin(t1188);
t1163 = cos(t1188);
t1252 = (t1110 * t1207 + t1209) * t1159 - (-t1106 * t1207 + t1210) * t1163;
t1011 = (t1067 * t1208 + t1068 * t1252) * t1054;
t1010 = t1011 ^ 2;
t1460 = t1154 * t1174 + t1184 * t1274;
t1485 = 0.2e1 * t1010 * ((t1114 * t1184 + t1157) * t1186 + t1309 + t1460);
t1227 = 0.1e1 / pkin(3);
t1371 = t1057 * t1227;
t1484 = t1371 * t1488;
t1374 = t1056 * t1227;
t1483 = t1374 * t1487;
t1377 = t1055 * t1227;
t1482 = t1377 * t1486;
t1380 = t1054 * t1227;
t1481 = t1380 * t1485;
t1257 = rSges(3,1) * t1202 - rSges(3,2) * t1196;
t1476 = m(3) * t1257;
t1258 = rSges(3,1) * t1200 - rSges(3,2) * t1194;
t1475 = m(3) * t1258;
t1259 = rSges(3,1) * t1198 - rSges(3,2) * t1192;
t1474 = m(3) * t1259;
t1260 = rSges(3,1) * t1186 - rSges(3,2) * t1184;
t1473 = m(3) * t1260;
t1116 = pkin(2) * t1187 + t1185 * t1211;
t1328 = t1183 * t1187;
t1331 = t1182 * t1183;
t1426 = pkin(3) * t1186;
t1042 = (t1180 * t1185 - t1182 * t1328) * t1426 - t1116 * t1331 + t1115 * t1180;
t1348 = t1180 * t1183;
t1043 = (t1180 * t1328 + t1182 * t1185) * t1426 + t1116 * t1348 + t1182 * t1115;
t1003 = (-t1042 * t1252 + t1043 * t1208) * t1380;
t1128 = rSges(3,1) * t1184 + rSges(3,2) * t1186;
t1204 = pkin(6) + rSges(3,3);
t1440 = m(3) * t1204;
t1137 = m(2) * rSges(2,2) - t1440;
t1150 = m(2) * rSges(2,1) + t1212;
t1175 = m(1) + m(2) + m(3);
t1370 = ((t1150 + t1473) * t1187 - t1185 * t1137) * t1181;
t1409 = t1011 * t1187;
t1415 = 0.2e1 * m(3);
t1441 = m(3) * t1183;
t1066 = -t1128 * t1346 + t1183 * t1260;
t1445 = m(3) * t1066;
t1434 = pkin(3) * t1003;
t1302 = t1184 * t1434;
t1317 = t1185 * t1186;
t1344 = t1181 * t1187;
t1347 = t1181 * t1184;
t1408 = t1011 * t1211;
t1414 = t1003 * t1054;
t1290 = t1184 * t1408;
t954 = t1290 - t1434;
t918 = ((t1003 * t1183 + t1011 * t1344) * t1430 + ((-t1302 + t1408) * t1185 + pkin(2) * t1409) * t1345 + t1183 * t954) * t1054 * t1011 + pkin(3) * (t1003 * t1344 + (t1174 * t1183 - t1317 * t1347 - t1183) * t1011) * t1414;
t1094 = pkin(3) * t1317 + t1115;
t1318 = t1183 * t1227;
t1228 = pkin(2) ^ 2;
t1155 = t1211 ^ 2 + t1228;
t1226 = pkin(3) ^ 2;
t1452 = 0.2e1 * pkin(2);
t1421 = pkin(3) * t1452;
t1420 = t1011 * (-t1211 * t1302 + (t1174 * t1226 + t1186 * t1421 + t1155) * t1011);
t922 = t1054 * t1318 * t1420 + (-t1183 * t1290 + (-t1094 * t1347 + t1183 * (pkin(2) * t1186 + t1430)) * t1003) / (t1094 * t1345 + (pkin(2) + t1426) * t1330) * t1003;
t930 = (-t1186 * t1420 - (pkin(2) * t1003 - t1186 * t954) * t1434) * t1054;
t999 = t1003 ^ 2;
t1472 = -t1175 * t930 - t1370 * t918 - t1445 * t922 + ((-t1010 * t1150 - (t1010 + t999) * t1473) * t1185 - (t1003 * t1128 * t1415 + t1011 * t1137) * t1409) * t1181 - t999 * t1128 * t1441;
t1120 = pkin(2) * t1199 + t1193 * t1211;
t1321 = t1183 * t1199;
t1425 = pkin(3) * t1198;
t1044 = (t1180 * t1193 - t1182 * t1321) * t1425 - t1120 * t1331 + t1117 * t1180;
t1047 = (t1180 * t1321 + t1182 * t1193) * t1425 + t1120 * t1348 + t1182 * t1117;
t1007 = (-t1044 * t1251 + t1047 * t1208) * t1377;
t1004 = t1007 ^ 2;
t1132 = rSges(3,1) * t1192 + rSges(3,2) * t1198;
t1369 = ((t1150 + t1474) * t1199 - t1193 * t1137) * t1181;
t1404 = t1015 * t1199;
t1075 = -t1132 * t1342 + t1183 * t1259;
t1444 = m(3) * t1075;
t1433 = pkin(3) * t1007;
t1301 = t1192 * t1433;
t1316 = t1193 * t1198;
t1336 = t1181 * t1199;
t1343 = t1181 * t1192;
t1403 = t1015 * t1211;
t1413 = t1007 * t1055;
t1289 = t1192 * t1403;
t956 = t1289 - t1433;
t919 = ((t1007 * t1183 + t1015 * t1336) * t1429 + ((-t1301 + t1403) * t1193 + pkin(2) * t1404) * t1337 + t1183 * t956) * t1055 * t1015 + pkin(3) * (t1007 * t1336 + (t1176 * t1183 - t1316 * t1343 - t1183) * t1015) * t1413;
t1095 = pkin(3) * t1316 + t1117;
t1419 = t1015 * (-t1211 * t1301 + (t1176 * t1226 + t1198 * t1421 + t1155) * t1015);
t923 = t1055 * t1318 * t1419 + (-t1183 * t1289 + (-t1095 * t1343 + (pkin(2) * t1198 + t1429) * t1183) * t1007) / (t1095 * t1337 + (pkin(2) + t1425) * t1327) * t1007;
t931 = (-t1198 * t1419 - (pkin(2) * t1007 - t1198 * t956) * t1433) * t1055;
t1471 = -t1175 * t931 - t1369 * t919 - t1444 * t923 + ((-t1012 * t1150 - (t1012 + t1004) * t1474) * t1193 - (t1007 * t1132 * t1415 + t1015 * t1137) * t1404) * t1181 - t1004 * t1132 * t1441;
t1121 = pkin(2) * t1201 + t1195 * t1211;
t1320 = t1183 * t1201;
t1424 = pkin(3) * t1200;
t1045 = (t1180 * t1195 - t1182 * t1320) * t1424 - t1121 * t1331 + t1118 * t1180;
t1048 = (t1180 * t1320 + t1182 * t1195) * t1424 + t1121 * t1348 + t1182 * t1118;
t1008 = (-t1045 * t1250 + t1048 * t1208) * t1374;
t1005 = t1008 ^ 2;
t1133 = rSges(3,1) * t1194 + rSges(3,2) * t1200;
t1368 = ((t1150 + t1475) * t1201 - t1195 * t1137) * t1181;
t1402 = t1016 * t1201;
t1076 = -t1133 * t1340 + t1183 * t1258;
t1443 = m(3) * t1076;
t1432 = pkin(3) * t1008;
t1300 = t1194 * t1432;
t1315 = t1195 * t1200;
t1334 = t1181 * t1201;
t1341 = t1181 * t1194;
t1401 = t1016 * t1211;
t1412 = t1008 * t1056;
t1288 = t1194 * t1401;
t957 = t1288 - t1432;
t920 = ((t1008 * t1183 + t1016 * t1334) * t1428 + ((-t1300 + t1401) * t1195 + pkin(2) * t1402) * t1335 + t1183 * t957) * t1056 * t1016 + pkin(3) * (t1008 * t1334 + (t1177 * t1183 - t1315 * t1341 - t1183) * t1016) * t1412;
t1096 = pkin(3) * t1315 + t1118;
t1418 = t1016 * (-t1211 * t1300 + (t1177 * t1226 + t1200 * t1421 + t1155) * t1016);
t924 = t1056 * t1318 * t1418 + (-t1183 * t1288 + (-t1096 * t1341 + (pkin(2) * t1200 + t1428) * t1183) * t1008) / (t1096 * t1335 + (pkin(2) + t1424) * t1325) * t1008;
t932 = (-t1200 * t1418 - (pkin(2) * t1008 - t1200 * t957) * t1432) * t1056;
t1470 = -t1175 * t932 - t1368 * t920 - t1443 * t924 + ((-t1013 * t1150 - (t1013 + t1005) * t1475) * t1195 - (t1008 * t1133 * t1415 + t1016 * t1137) * t1402) * t1181 - t1005 * t1133 * t1441;
t1122 = pkin(2) * t1203 + t1197 * t1211;
t1319 = t1183 * t1203;
t1423 = pkin(3) * t1202;
t1046 = (t1180 * t1197 - t1182 * t1319) * t1423 - t1122 * t1331 + t1119 * t1180;
t1049 = (t1180 * t1319 + t1182 * t1197) * t1423 + t1122 * t1348 + t1182 * t1119;
t1009 = (-t1046 * t1249 + t1049 * t1208) * t1371;
t1006 = t1009 ^ 2;
t1134 = rSges(3,1) * t1196 + rSges(3,2) * t1202;
t1367 = ((t1150 + t1476) * t1203 - t1197 * t1137) * t1181;
t1400 = t1017 * t1203;
t1077 = -t1134 * t1338 + t1183 * t1257;
t1442 = m(3) * t1077;
t1431 = pkin(3) * t1009;
t1299 = t1196 * t1431;
t1314 = t1197 * t1202;
t1332 = t1181 * t1203;
t1339 = t1181 * t1196;
t1399 = t1017 * t1211;
t1411 = t1009 * t1057;
t1287 = t1196 * t1399;
t958 = t1287 - t1431;
t921 = ((t1009 * t1183 + t1017 * t1332) * t1427 + ((-t1299 + t1399) * t1197 + pkin(2) * t1400) * t1333 + t1183 * t958) * t1057 * t1017 + pkin(3) * (t1009 * t1332 + (t1178 * t1183 - t1314 * t1339 - t1183) * t1017) * t1411;
t1097 = pkin(3) * t1314 + t1119;
t1417 = t1017 * (-t1211 * t1299 + (t1178 * t1226 + t1202 * t1421 + t1155) * t1017);
t925 = t1057 * t1318 * t1417 + (-t1183 * t1287 + (-t1097 * t1339 + (pkin(2) * t1202 + t1427) * t1183) * t1009) / (t1097 * t1333 + (pkin(2) + t1423) * t1323) * t1009;
t933 = (-t1202 * t1417 - (pkin(2) * t1009 - t1202 * t958) * t1431) * t1057;
t1469 = -t1175 * t933 - t1367 * t921 - t1442 * t925 + ((-t1014 * t1150 - (t1014 + t1006) * t1476) * t1197 - (t1009 * t1134 * t1415 + t1017 * t1137) * t1400) * t1181 - t1006 * t1134 * t1441;
t1464 = t1054 * t1472;
t1463 = t1055 * t1471;
t1462 = t1056 * t1470;
t1461 = t1057 * t1469;
t1277 = t1046 * t1371;
t1446 = -t1154 / 0.2e1;
t1142 = rSges(3,1) * t1440 - Icges(3,5);
t1447 = t1142 / 0.4e1;
t1141 = rSges(3,2) * t1440 - Icges(3,6);
t1448 = -t1141 / 0.4e1;
t1135 = (-t1224 + t1225) * m(3) + Icges(3,2) - Icges(3,1);
t1449 = t1135 / 0.2e1;
t937 = (t1196 * t1448 + t1202 * t1447) * t1009 + ((t1196 * t1449 + t1157) * t1202 + t1446 + t1457) * t1017;
t1295 = t937 * t1411;
t1373 = t1057 * t1074;
t1138 = t1204 ^ 2 + t1224 + t1228;
t1168 = t1450 * t1452;
t1265 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) + Icges(2,3);
t1416 = -0.2e1 * rSges(3,2) * pkin(2);
t1451 = -0.2e1 * t1154;
t1041 = t1135 * t1178 + (t1196 * t1451 + t1168) * t1202 + (t1196 * t1416 + t1138) * m(3) + t1265;
t1085 = -t1141 * t1202 - t1196 * t1142;
t911 = -t1041 * t921 - t1085 * t925 - t1367 * t933;
t1136 = (t1224 + t1225) * m(3) + Icges(3,3);
t917 = -t1085 * t921 - t1136 * t925 - t1442 * t933;
t1456 = t1046 * t1484 + 0.4e1 * t1074 * t1295 + t917 * t1277 - t911 * t1373;
t1280 = t1045 * t1374;
t936 = (t1194 * t1448 + t1200 * t1447) * t1008 + ((t1194 * t1449 + t1157) * t1200 + t1446 + t1458) * t1016;
t1296 = t936 * t1412;
t1376 = t1056 * t1073;
t1040 = t1135 * t1177 + (t1194 * t1451 + t1168) * t1200 + (t1194 * t1416 + t1138) * m(3) + t1265;
t1084 = -t1141 * t1200 - t1194 * t1142;
t910 = -t1040 * t920 - t1084 * t924 - t1368 * t932;
t916 = -t1084 * t920 - t1136 * t924 - t1443 * t932;
t1455 = t1045 * t1483 + 0.4e1 * t1073 * t1296 + t916 * t1280 - t910 * t1376;
t1286 = t1042 * t1380;
t934 = (t1184 * t1448 + t1186 * t1447) * t1003 + ((t1184 * t1449 + t1157) * t1186 + t1446 + t1460) * t1011;
t1298 = t934 * t1414;
t1382 = t1054 * t1068;
t1038 = t1135 * t1174 + (t1184 * t1451 + t1168) * t1186 + (t1184 * t1416 + t1138) * m(3) + t1265;
t1082 = -t1141 * t1186 - t1184 * t1142;
t906 = -t1038 * t918 - t1082 * t922 - t1370 * t930;
t908 = -t1082 * t918 - t1136 * t922 - t1445 * t930;
t1454 = t1042 * t1481 + 0.4e1 * t1068 * t1298 + t908 * t1286 - t906 * t1382;
t1283 = t1044 * t1377;
t935 = (t1192 * t1448 + t1198 * t1447) * t1007 + ((t1192 * t1449 + t1157) * t1198 + t1446 + t1459) * t1015;
t1297 = t935 * t1413;
t1379 = t1055 * t1072;
t1039 = t1135 * t1176 + (t1192 * t1451 + t1168) * t1198 + (t1192 * t1416 + t1138) * m(3) + t1265;
t1083 = -t1141 * t1198 - t1192 * t1142;
t909 = -t1039 * t919 - t1083 * t923 - t1369 * t931;
t915 = -t1083 * t919 - t1136 * t923 - t1444 * t931;
t1453 = t1044 * t1482 + 0.4e1 * t1072 * t1297 + t915 * t1283 - t909 * t1379;
t1179 = t1207 ^ 2;
t1439 = m(3) * t1227;
t1438 = pkin(2) * t1184;
t1437 = pkin(2) * t1192;
t1436 = pkin(2) * t1194;
t1435 = pkin(2) * t1196;
t1422 = m(4) * t1179;
t1386 = t1042 * t1227;
t1385 = t1044 * t1227;
t1384 = t1045 * t1227;
t1383 = t1046 * t1227;
t1381 = t1054 * t1179;
t1378 = t1055 * t1179;
t1375 = t1056 * t1179;
t1372 = t1057 * t1179;
t1366 = t1068 * t1159;
t1365 = t1068 * t1163;
t1364 = t1072 * t1160;
t1363 = t1072 * t1164;
t1362 = t1073 * t1161;
t1361 = t1073 * t1165;
t1360 = t1074 * t1162;
t1359 = t1074 * t1166;
t1358 = t1082 * t1227;
t1357 = t1083 * t1227;
t1356 = t1084 * t1227;
t1355 = t1085 * t1227;
t1354 = t1136 * t1227;
t1349 = t1180 * t1181;
t1306 = t1066 * t1439;
t1305 = t1075 * t1439;
t1304 = t1076 * t1439;
t1303 = t1077 * t1439;
t1285 = t1159 * t1386;
t1284 = t1163 * t1386;
t1282 = t1160 * t1385;
t1281 = t1164 * t1385;
t1279 = t1161 * t1384;
t1278 = t1165 * t1384;
t1276 = t1162 * t1383;
t1275 = t1166 * t1383;
t1248 = t1106 * t1163 + t1110 * t1159;
t1247 = t1107 * t1164 + t1111 * t1160;
t1246 = t1108 * t1165 + t1112 * t1161;
t1245 = t1109 * t1166 + t1113 * t1162;
t1244 = pkin(3) * t1347 - t1115 * t1183;
t1243 = pkin(3) * t1343 - t1117 * t1183;
t1242 = pkin(3) * t1341 - t1118 * t1183;
t1241 = pkin(3) * t1339 - t1119 * t1183;
t1240 = -t1038 * t1068 + t1042 * t1358;
t1239 = -t1039 * t1072 + t1044 * t1357;
t1238 = -t1040 * t1073 + t1045 * t1356;
t1237 = -t1041 * t1074 + t1046 * t1355;
t1236 = t1042 * t1354 - t1068 * t1082;
t1235 = t1044 * t1354 - t1072 * t1083;
t1234 = t1045 * t1354 - t1073 * t1084;
t1233 = t1046 * t1354 - t1074 * t1085;
t1232 = t1042 * t1306 - t1068 * t1370;
t1231 = t1044 * t1305 - t1072 * t1369;
t1230 = t1045 * t1304 - t1073 * t1368;
t1229 = t1046 * t1303 - t1074 * t1367;
t1215 = rSges(4,1);
t1214 = rSges(4,2);
t1061 = t1122 * t1182 + t1180 * t1241;
t1060 = t1121 * t1182 + t1180 * t1242;
t1059 = t1120 * t1182 + t1180 * t1243;
t1058 = t1116 * t1182 + t1180 * t1244;
t1037 = -t1105 * t1427 - t1122 * t1180 * t1202 + (pkin(2) * t1339 + t1202 * t1241) * t1182;
t1036 = -t1104 * t1428 - t1121 * t1180 * t1200 + (pkin(2) * t1341 + t1200 * t1242) * t1182;
t1035 = -t1103 * t1429 - t1120 * t1180 * t1198 + (pkin(2) * t1343 + t1198 * t1243) * t1182;
t1034 = -t1099 * t1430 - t1116 * t1180 * t1186 + (pkin(2) * t1347 + t1186 * t1244) * t1182;
t1033 = -(t1102 * t1166 - t1162 * t1338) * t1427 + (t1061 * t1166 + t1081 * t1162) * t1202 + (t1162 * t1183 + t1166 * t1349) * t1435;
t1032 = (t1102 * t1162 + t1166 * t1338) * t1427 + (-t1061 * t1162 + t1081 * t1166) * t1202 + (-t1162 * t1349 + t1166 * t1183) * t1435;
t1031 = -(t1101 * t1165 - t1161 * t1340) * t1428 + (t1060 * t1165 + t1080 * t1161) * t1200 + (t1161 * t1183 + t1165 * t1349) * t1436;
t1030 = (t1101 * t1161 + t1165 * t1340) * t1428 + (-t1060 * t1161 + t1080 * t1165) * t1200 + (-t1161 * t1349 + t1165 * t1183) * t1436;
t1029 = -(t1100 * t1164 - t1160 * t1342) * t1429 + (t1059 * t1164 + t1079 * t1160) * t1198 + (t1160 * t1183 + t1164 * t1349) * t1437;
t1028 = (t1100 * t1160 + t1164 * t1342) * t1429 + (-t1059 * t1160 + t1079 * t1164) * t1198 + (-t1160 * t1349 + t1164 * t1183) * t1437;
t1027 = -(t1098 * t1163 - t1159 * t1346) * t1430 + (t1058 * t1163 + t1078 * t1159) * t1186 + (t1159 * t1183 + t1163 * t1349) * t1438;
t1026 = (t1098 * t1159 + t1163 * t1346) * t1430 + (-t1058 * t1159 + t1078 * t1163) * t1186 + (-t1159 * t1349 + t1163 * t1183) * t1438;
t1025 = t1245 * t1373;
t1024 = t1246 * t1376;
t1023 = t1247 * t1379;
t1022 = t1248 * t1382;
t1021 = t1245 * t1277;
t1020 = t1246 * t1280;
t1019 = t1247 * t1283;
t1018 = t1248 * t1286;
t1002 = (t1037 * t1442 + t1049 * t1354 + t1071 * t1085) * t1057;
t1001 = (t1036 * t1443 + t1048 * t1354 + t1070 * t1084) * t1056;
t1000 = (t1035 * t1444 + t1047 * t1354 + t1069 * t1083) * t1055;
t998 = (t1037 * t1175 + t1049 * t1303 + t1071 * t1367) * t1057;
t997 = (t1036 * t1175 + t1048 * t1304 + t1070 * t1368) * t1056;
t996 = (t1035 * t1175 + t1047 * t1305 + t1069 * t1369) * t1055;
t995 = (t1034 * t1445 + t1043 * t1354 + t1067 * t1082) * t1054;
t994 = (t1034 * t1175 + t1043 * t1306 + t1067 * t1370) * t1054;
t993 = (t1032 * t1113 - t1033 * t1109) * t1057;
t992 = (t1030 * t1112 - t1031 * t1108) * t1056;
t991 = (t1028 * t1111 - t1029 * t1107) * t1055;
t990 = (t1026 * t1110 - t1027 * t1106) * t1054;
t989 = (t1037 * t1367 + t1041 * t1071 + t1049 * t1355) * t1057;
t988 = (t1036 * t1368 + t1040 * t1070 + t1048 * t1356) * t1056;
t987 = (t1035 * t1369 + t1039 * t1069 + t1047 * t1357) * t1055;
t986 = (t1034 * t1370 + t1038 * t1067 + t1043 * t1358) * t1054;
t985 = (t1033 * t1442 + t1166 * t1233) * t1057;
t984 = (t1032 * t1442 - t1162 * t1233) * t1057;
t983 = (t1031 * t1443 + t1165 * t1234) * t1056;
t982 = (t1030 * t1443 - t1161 * t1234) * t1056;
t981 = (t1029 * t1444 + t1164 * t1235) * t1055;
t980 = (t1028 * t1444 - t1160 * t1235) * t1055;
t979 = (t1033 * t1175 + t1166 * t1229) * t1057;
t978 = (t1032 * t1175 - t1162 * t1229) * t1057;
t977 = (t1031 * t1175 + t1165 * t1230) * t1056;
t976 = (t1030 * t1175 - t1161 * t1230) * t1056;
t975 = (t1029 * t1175 + t1164 * t1231) * t1055;
t974 = (t1028 * t1175 - t1160 * t1231) * t1055;
t973 = (t1027 * t1445 + t1163 * t1236) * t1054;
t972 = (t1026 * t1445 - t1159 * t1236) * t1054;
t971 = (t1027 * t1175 + t1163 * t1232) * t1054;
t970 = (t1026 * t1175 - t1159 * t1232) * t1054;
t969 = (t1033 * t1367 + t1166 * t1237) * t1057;
t968 = (t1032 * t1367 - t1162 * t1237) * t1057;
t967 = (t1031 * t1368 + t1165 * t1238) * t1056;
t966 = (t1030 * t1368 - t1161 * t1238) * t1056;
t965 = (t1029 * t1369 + t1164 * t1239) * t1055;
t964 = (t1028 * t1369 - t1160 * t1239) * t1055;
t963 = (t1027 * t1370 + t1163 * t1240) * t1054;
t962 = (t1026 * t1370 - t1159 * t1240) * t1054;
t953 = -t1021 * t1136 + t1025 * t1085 + t1442 * t993;
t952 = -t1020 * t1136 + t1024 * t1084 + t1443 * t992;
t951 = -t1019 * t1136 + t1023 * t1083 + t1444 * t991;
t950 = -t1021 * t1442 + t1025 * t1367 + t1175 * t993;
t949 = -t1020 * t1443 + t1024 * t1368 + t1175 * t992;
t948 = -t1019 * t1444 + t1023 * t1369 + t1175 * t991;
t947 = -t1018 * t1136 + t1022 * t1082 + t1445 * t990;
t946 = -t1018 * t1445 + t1022 * t1370 + t1175 * t990;
t945 = -t1021 * t1085 + t1025 * t1041 + t1367 * t993;
t944 = -t1020 * t1084 + t1024 * t1040 + t1368 * t992;
t943 = -t1019 * t1083 + t1023 * t1039 + t1369 * t991;
t942 = -t1018 * t1082 + t1022 * t1038 + t1370 * t990;
t1 = [(-(t1033 * t979 + t1275 * t985 - t1359 * t969) * t1113 - (t1032 * t979 - t1276 * t985 + t1360 * t969) * t1109) * t1372 + (-(t1031 * t977 + t1278 * t983 - t1361 * t967) * t1112 - (t1030 * t977 - t1279 * t983 + t1362 * t967) * t1108) * t1375 + (-(t1029 * t975 + t1281 * t981 - t1363 * t965) * t1111 - (t1028 * t975 - t1282 * t981 + t1364 * t965) * t1107) * t1378 + (t1172 * t1214 - t1173 * t1215) * t1422 + (-(t1027 * t971 + t1284 * t973 - t1365 * t963) * t1110 - (t1026 * t971 - t1285 * t973 + t1366 * t963) * t1106) * t1381 + t1027 * t1464 + t1029 * t1463 + t1031 * t1462 + t1033 * t1461 + t1456 * t1166 + t1455 * t1165 + t1453 * t1164 + t1454 * t1163; (-(t1033 * t978 + t1275 * t984 - t1359 * t968) * t1113 - (t1032 * t978 - t1276 * t984 + t1360 * t968) * t1109) * t1372 + (-(t1031 * t976 + t1278 * t982 - t1361 * t966) * t1112 - (t1030 * t976 - t1279 * t982 + t1362 * t966) * t1108) * t1375 + (-(t1029 * t974 + t1281 * t980 - t1363 * t964) * t1111 - (t1028 * t974 - t1282 * t980 + t1364 * t964) * t1107) * t1378 + (-(t1027 * t970 + t1284 * t972 - t1365 * t962) * t1110 - (t1026 * t970 - t1285 * t972 + t1366 * t962) * t1106) * t1381 - (t1172 * t1215 + t1173 * t1214) * t1422 + t1026 * t1464 + t1028 * t1463 + t1030 * t1462 + t1032 * t1461 - t1456 * t1162 - t1455 * t1161 - t1453 * t1160 - t1454 * t1159; (-(t1002 * t1275 + t1033 * t998 - t1359 * t989) * t1113 - (-t1002 * t1276 + t1032 * t998 + t1360 * t989) * t1109) * t1372 + (-(t1001 * t1278 + t1031 * t997 - t1361 * t988) * t1112 - (-t1001 * t1279 + t1030 * t997 + t1362 * t988) * t1108) * t1375 + (-(t1000 * t1281 + t1029 * t996 - t1363 * t987) * t1111 - (-t1000 * t1282 + t1028 * t996 + t1364 * t987) * t1107) * t1378 + (-(t1027 * t994 + t1284 * t995 - t1365 * t986) * t1110 - (t1026 * t994 - t1285 * t995 + t1366 * t986) * t1106) * t1381 + t1034 * t1464 + t1035 * t1463 + t1036 * t1462 + t1037 * t1461 + (t1057 * t911 - 0.4e1 * t1295) * t1071 + (t1056 * t910 - 0.4e1 * t1296) * t1070 + (t1055 * t909 - 0.4e1 * t1297) * t1069 + (t1054 * t906 - 0.4e1 * t1298) * t1067 + (t1371 * t917 + t1484) * t1049 + (t1374 * t916 + t1483) * t1048 + (t1377 * t915 + t1482) * t1047 + (t1380 * t908 + t1481) * t1043; (-(t1033 * t950 + t1275 * t953 - t1359 * t945) * t1113 - (t1032 * t950 - t1276 * t953 + t1360 * t945) * t1109) * t1372 + (-(t1031 * t949 + t1278 * t952 - t1361 * t944) * t1112 - (t1030 * t949 - t1279 * t952 + t1362 * t944) * t1108) * t1375 + (-(t1029 * t948 + t1281 * t951 - t1363 * t943) * t1111 - (t1028 * t948 - t1282 * t951 + t1364 * t943) * t1107) * t1378 + (-(t1027 * t946 + t1284 * t947 - t1365 * t942) * t1110 - (t1026 * t946 - t1285 * t947 + t1366 * t942) * t1106) * t1381 + t1469 * t993 + t1470 * t992 + t1471 * t991 + t1472 * t990 + (-0.4e1 * t1009 * t937 + t911) * t1025 + (-0.4e1 * t1008 * t936 + t910) * t1024 + (-0.4e1 * t1007 * t935 + t909) * t1023 + (-0.4e1 * t1003 * t934 + t906) * t1022 - (t1488 + t917) * t1021 - (t1487 + t916) * t1020 - (t1486 + t915) * t1019 - (t1485 + t908) * t1018;];
taucX  = t1;
