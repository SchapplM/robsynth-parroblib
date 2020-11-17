% Calculate vector of centrifugal and coriolis load on the joints for
% P4PRRRR8V2G1A0
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
% Datum: 2020-08-07 11:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:09:35
% EndTime: 2020-08-07 11:09:48
% DurationCPUTime: 13.67s
% Computational Cost: add. (57452->595), mult. (109766->1089), div. (5368->13), fcn. (110828->30), ass. (0->431)
t1198 = legFrame(4,3);
t1172 = sin(t1198);
t1176 = cos(t1198);
t1194 = sin(pkin(8));
t1196 = cos(pkin(8));
t1088 = -t1172 * t1194 + t1176 * t1196;
t1092 = t1172 * t1196 + t1176 * t1194;
t1204 = cos(qJ(3,4));
t1159 = pkin(3) * t1204 + pkin(2);
t1205 = cos(qJ(2,4));
t1224 = pkin(7) + pkin(6);
t1160 = t1205 * t1224;
t1203 = sin(qJ(2,4));
t1105 = t1159 * t1203 - t1160;
t1197 = cos(pkin(4));
t1310 = t1224 * t1203;
t1352 = (t1159 * t1205 + t1310) * t1197;
t1036 = -t1088 * t1352 + t1092 * t1105;
t1037 = -t1088 * t1105 - t1092 * t1352;
t1226 = xP(4);
t1186 = sin(t1226);
t1187 = cos(t1226);
t1229 = koppelP(4,2);
t1233 = koppelP(4,1);
t1116 = -t1186 * t1229 + t1187 * t1233;
t1221 = xDP(4);
t1222 = xDP(2);
t1076 = t1116 * t1221 + t1222;
t1112 = t1186 * t1233 + t1187 * t1229;
t1223 = xDP(1);
t1080 = -t1112 * t1221 + t1223;
t1202 = sin(qJ(3,4));
t1323 = t1197 * t1202;
t1120 = t1159 * t1323;
t1195 = sin(pkin(4));
t1337 = t1195 * t1204;
t1060 = 0.1e1 / (t1105 * t1337 + t1120);
t1240 = 0.1e1 / pkin(3);
t1360 = t1060 * t1240;
t1009 = (t1036 * t1080 + t1037 * t1076) * t1360;
t1237 = rSges(3,2) ^ 2;
t1238 = rSges(3,1) ^ 2;
t1145 = (-t1237 + t1238) * m(3) + Icges(3,2) - Icges(3,1);
t1218 = pkin(6) + rSges(3,3);
t1241 = pkin(2) ^ 2;
t1148 = t1218 ^ 2 + t1237 + t1241;
t1415 = m(3) * rSges(3,1);
t1417 = 0.2e1 * pkin(2);
t1182 = t1415 * t1417;
t1188 = t1204 ^ 2;
t1258 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) + Icges(2,3);
t1380 = -0.2e1 * rSges(3,2) * pkin(2);
t1301 = rSges(3,2) * t1415;
t1168 = -Icges(3,4) + t1301;
t1416 = -0.2e1 * t1168;
t1044 = t1145 * t1188 + (t1202 * t1416 + t1182) * t1204 + (t1202 * t1380 + t1148) * m(3) + t1258;
t1405 = m(3) * t1218;
t1151 = rSges(3,2) * t1405 - Icges(3,6);
t1152 = rSges(3,1) * t1405 - Icges(3,5);
t1072 = -t1151 * t1204 - t1152 * t1202;
t1225 = pkin(2) * m(3);
t1300 = t1225 / 0.2e1;
t1171 = rSges(3,2) * t1300;
t1147 = m(2) * rSges(2,2) - t1405;
t1161 = m(2) * rSges(2,1) + t1225;
t1253 = rSges(3,1) * t1204 - rSges(3,2) * t1202;
t1446 = t1253 * m(3);
t1064 = (t1161 + t1446) * t1205 - t1147 * t1203;
t1356 = t1064 * t1195;
t1411 = -t1168 / 0.2e1;
t1412 = t1152 / 0.4e1;
t1413 = -t1151 / 0.4e1;
t1414 = t1145 / 0.2e1;
t1267 = rSges(3,1) * t1300;
t1421 = t1168 * t1188 + t1202 * t1267;
t1125 = pkin(2) * t1203 - t1160;
t1338 = t1195 * t1203;
t1400 = pkin(3) * t1188;
t1052 = 0.1e1 / ((pkin(3) * t1323 + t1125 * t1195) * t1204 + pkin(2) * t1323 + t1338 * t1400);
t1404 = pkin(3) * t1009;
t1287 = t1202 * t1404;
t1314 = t1203 * t1204;
t1336 = t1195 * t1205;
t1339 = t1195 * t1202;
t1322 = t1197 * t1203;
t1028 = -t1088 * t1337 - (t1088 * t1322 + t1092 * t1205) * t1202;
t1029 = -t1092 * t1337 - (-t1088 * t1205 + t1092 * t1322) * t1202;
t997 = (t1028 * t1080 + t1029 * t1076) * t1052;
t1381 = t1224 * t997;
t1382 = t1205 * t997;
t1279 = t1202 * t1381;
t968 = t1279 - t1404;
t916 = (((t1009 * t1197 + t1336 * t997) * t1400 + ((-t1287 + t1381) * t1203 + pkin(2) * t1382) * t1337 + t968 * t1197) * t997 + (t1009 * t1336 + (t1188 * t1197 - t1314 * t1339 - t1197) * t997) * t1404) * t1052;
t1084 = pkin(3) * t1314 + t1125;
t1315 = t1197 * t1240;
t1169 = t1224 ^ 2 + t1241;
t1239 = pkin(3) ^ 2;
t1394 = pkin(3) * t1417;
t1395 = (-t1224 * t1287 + (t1188 * t1239 + t1204 * t1394 + t1169) * t997) * t997;
t920 = t1052 * t1315 * t1395 + (-t1197 * t1279 + (-t1084 * t1339 + (pkin(2) * t1204 + t1400) * t1197) * t1009) / (t1084 * t1337 + t1120) * t1009;
t929 = (-t1204 * t1395 - (pkin(2) * t1009 - t1204 * t968) * t1404) * t1052;
t1461 = -t1044 * t916 - t1072 * t920 - t1356 * t929 - 0.4e1 * t1009 * ((t1413 * t1202 + t1412 * t1204) * t1009 + ((t1202 * t1414 + t1171) * t1204 + t1411 + t1421) * t997);
t1124 = (t1238 / 0.2e1 - t1237 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t1146 = (t1237 + t1238) * m(3) + Icges(3,3);
t1302 = -t1301 / 0.2e1 + Icges(3,4) / 0.2e1;
t1138 = rSges(3,1) * t1202 + rSges(3,2) * t1204;
t1068 = -t1138 * t1338 + t1253 * t1197;
t1410 = m(3) * t1068;
t996 = t997 ^ 2;
t1460 = -t1072 * t916 - t1146 * t920 - t1410 * t929 + 0.2e1 * ((t1124 * t1202 + t1171) * t1204 + t1302 + t1421) * t996;
t1199 = legFrame(3,3);
t1173 = sin(t1199);
t1177 = cos(t1199);
t1089 = -t1173 * t1194 + t1177 * t1196;
t1093 = t1173 * t1196 + t1177 * t1194;
t1206 = sin(qJ(3,3));
t1213 = cos(qJ(2,3));
t1207 = sin(qJ(2,3));
t1320 = t1197 * t1207;
t1212 = cos(qJ(3,3));
t1329 = t1195 * t1212;
t1030 = -t1089 * t1329 - (t1089 * t1320 + t1093 * t1213) * t1206;
t1033 = -t1093 * t1329 - (-t1089 * t1213 + t1093 * t1320) * t1206;
t1165 = t1213 * t1224;
t1127 = pkin(2) * t1207 - t1165;
t1321 = t1197 * t1206;
t1334 = t1195 * t1207;
t1190 = t1212 ^ 2;
t1399 = pkin(3) * t1190;
t1053 = 0.1e1 / ((pkin(3) * t1321 + t1127 * t1195) * t1212 + pkin(2) * t1321 + t1334 * t1399);
t1230 = koppelP(3,2);
t1234 = koppelP(3,1);
t1117 = -t1186 * t1230 + t1187 * t1234;
t1077 = t1117 * t1221 + t1222;
t1113 = t1186 * t1234 + t1187 * t1230;
t1081 = -t1113 * t1221 + t1223;
t1001 = (t1030 * t1081 + t1033 * t1077) * t1053;
t1162 = pkin(3) * t1212 + pkin(2);
t1109 = t1162 * t1207 - t1165;
t1309 = t1224 * t1207;
t1351 = (t1162 * t1213 + t1309) * t1197;
t1038 = -t1089 * t1351 + t1093 * t1109;
t1041 = -t1089 * t1109 - t1093 * t1351;
t1121 = t1162 * t1321;
t1061 = 0.1e1 / (t1109 * t1329 + t1121);
t1359 = t1061 * t1240;
t1013 = (t1038 * t1081 + t1041 * t1077) * t1359;
t1045 = t1145 * t1190 + (t1206 * t1416 + t1182) * t1212 + (t1206 * t1380 + t1148) * m(3) + t1258;
t1073 = -t1151 * t1212 - t1152 * t1206;
t1252 = rSges(3,1) * t1212 - rSges(3,2) * t1206;
t1447 = t1252 * m(3);
t1065 = (t1161 + t1447) * t1213 - t1147 * t1207;
t1355 = t1065 * t1195;
t1420 = t1168 * t1190 + t1206 * t1267;
t1403 = pkin(3) * t1013;
t1286 = t1206 * t1403;
t1313 = t1207 * t1212;
t1328 = t1195 * t1213;
t1335 = t1195 * t1206;
t1377 = t1001 * t1224;
t1378 = t1001 * t1213;
t1278 = t1206 * t1377;
t970 = t1278 - t1403;
t917 = (((t1001 * t1328 + t1013 * t1197) * t1399 + ((-t1286 + t1377) * t1207 + pkin(2) * t1378) * t1329 + t970 * t1197) * t1001 + (t1013 * t1328 + (t1190 * t1197 - t1313 * t1335 - t1197) * t1001) * t1403) * t1053;
t1085 = pkin(3) * t1313 + t1127;
t1393 = t1001 * (-t1224 * t1286 + (t1190 * t1239 + t1212 * t1394 + t1169) * t1001);
t921 = t1053 * t1315 * t1393 + (-t1197 * t1278 + (-t1085 * t1335 + t1197 * (pkin(2) * t1212 + t1399)) * t1013) / (t1085 * t1329 + t1121) * t1013;
t933 = (-t1212 * t1393 - (pkin(2) * t1013 - t1212 * t970) * t1403) * t1053;
t1459 = -t1045 * t917 - t1073 * t921 - t1355 * t933 - 0.4e1 * t1013 * ((t1413 * t1206 + t1412 * t1212) * t1013 + ((t1206 * t1414 + t1171) * t1212 + t1411 + t1420) * t1001);
t1200 = legFrame(2,3);
t1174 = sin(t1200);
t1178 = cos(t1200);
t1090 = -t1174 * t1194 + t1178 * t1196;
t1094 = t1174 * t1196 + t1178 * t1194;
t1208 = sin(qJ(3,2));
t1215 = cos(qJ(2,2));
t1209 = sin(qJ(2,2));
t1318 = t1197 * t1209;
t1214 = cos(qJ(3,2));
t1327 = t1195 * t1214;
t1031 = -t1090 * t1327 - (t1090 * t1318 + t1094 * t1215) * t1208;
t1034 = -t1094 * t1327 - (-t1090 * t1215 + t1094 * t1318) * t1208;
t1166 = t1215 * t1224;
t1128 = pkin(2) * t1209 - t1166;
t1319 = t1197 * t1208;
t1332 = t1195 * t1209;
t1191 = t1214 ^ 2;
t1398 = pkin(3) * t1191;
t1054 = 0.1e1 / ((pkin(3) * t1319 + t1128 * t1195) * t1214 + pkin(2) * t1319 + t1332 * t1398);
t1231 = koppelP(2,2);
t1235 = koppelP(2,1);
t1118 = -t1186 * t1231 + t1187 * t1235;
t1078 = t1118 * t1221 + t1222;
t1114 = t1186 * t1235 + t1187 * t1231;
t1082 = -t1114 * t1221 + t1223;
t1002 = (t1031 * t1082 + t1034 * t1078) * t1054;
t1163 = pkin(3) * t1214 + pkin(2);
t1110 = t1163 * t1209 - t1166;
t1308 = t1224 * t1209;
t1350 = (t1163 * t1215 + t1308) * t1197;
t1039 = -t1090 * t1350 + t1094 * t1110;
t1042 = -t1090 * t1110 - t1094 * t1350;
t1122 = t1163 * t1319;
t1062 = 0.1e1 / (t1110 * t1327 + t1122);
t1358 = t1062 * t1240;
t1014 = (t1039 * t1082 + t1042 * t1078) * t1358;
t1046 = t1145 * t1191 + (t1208 * t1416 + t1182) * t1214 + (t1208 * t1380 + t1148) * m(3) + t1258;
t1074 = -t1151 * t1214 - t1152 * t1208;
t1251 = rSges(3,1) * t1214 - rSges(3,2) * t1208;
t1448 = t1251 * m(3);
t1066 = (t1161 + t1448) * t1215 - t1147 * t1209;
t1354 = t1066 * t1195;
t1419 = t1168 * t1191 + t1208 * t1267;
t1402 = pkin(3) * t1014;
t1285 = t1208 * t1402;
t1312 = t1209 * t1214;
t1326 = t1195 * t1215;
t1333 = t1195 * t1208;
t1375 = t1002 * t1224;
t1376 = t1002 * t1215;
t1277 = t1208 * t1375;
t971 = t1277 - t1402;
t918 = (((t1002 * t1326 + t1014 * t1197) * t1398 + ((-t1285 + t1375) * t1209 + pkin(2) * t1376) * t1327 + t971 * t1197) * t1002 + (t1014 * t1326 + (t1191 * t1197 - t1312 * t1333 - t1197) * t1002) * t1402) * t1054;
t1086 = pkin(3) * t1312 + t1128;
t1392 = t1002 * (-t1224 * t1285 + (t1191 * t1239 + t1214 * t1394 + t1169) * t1002);
t922 = t1054 * t1315 * t1392 + (-t1197 * t1277 + (-t1086 * t1333 + t1197 * (pkin(2) * t1214 + t1398)) * t1014) / (t1086 * t1327 + t1122) * t1014;
t934 = (-t1214 * t1392 - (pkin(2) * t1014 - t1214 * t971) * t1402) * t1054;
t1458 = -t1046 * t918 - t1074 * t922 - t1354 * t934 - 0.4e1 * t1014 * ((t1413 * t1208 + t1412 * t1214) * t1014 + ((t1208 * t1414 + t1171) * t1214 + t1411 + t1419) * t1002);
t1201 = legFrame(1,3);
t1175 = sin(t1201);
t1179 = cos(t1201);
t1091 = -t1175 * t1194 + t1179 * t1196;
t1095 = t1175 * t1196 + t1179 * t1194;
t1210 = sin(qJ(3,1));
t1217 = cos(qJ(2,1));
t1211 = sin(qJ(2,1));
t1316 = t1197 * t1211;
t1216 = cos(qJ(3,1));
t1325 = t1195 * t1216;
t1032 = -t1091 * t1325 - (t1091 * t1316 + t1095 * t1217) * t1210;
t1035 = -t1095 * t1325 - (-t1091 * t1217 + t1095 * t1316) * t1210;
t1167 = t1217 * t1224;
t1129 = pkin(2) * t1211 - t1167;
t1317 = t1197 * t1210;
t1330 = t1195 * t1211;
t1192 = t1216 ^ 2;
t1397 = pkin(3) * t1192;
t1055 = 0.1e1 / ((pkin(3) * t1317 + t1129 * t1195) * t1216 + pkin(2) * t1317 + t1330 * t1397);
t1232 = koppelP(1,2);
t1236 = koppelP(1,1);
t1119 = -t1186 * t1232 + t1187 * t1236;
t1079 = t1119 * t1221 + t1222;
t1115 = t1186 * t1236 + t1187 * t1232;
t1083 = -t1115 * t1221 + t1223;
t1003 = (t1032 * t1083 + t1035 * t1079) * t1055;
t1164 = pkin(3) * t1216 + pkin(2);
t1111 = t1164 * t1211 - t1167;
t1307 = t1224 * t1211;
t1349 = (t1164 * t1217 + t1307) * t1197;
t1040 = -t1091 * t1349 + t1095 * t1111;
t1043 = -t1091 * t1111 - t1095 * t1349;
t1123 = t1164 * t1317;
t1063 = 0.1e1 / (t1111 * t1325 + t1123);
t1357 = t1063 * t1240;
t1015 = (t1040 * t1083 + t1043 * t1079) * t1357;
t1047 = t1145 * t1192 + (t1210 * t1416 + t1182) * t1216 + (t1210 * t1380 + t1148) * m(3) + t1258;
t1075 = -t1151 * t1216 - t1152 * t1210;
t1250 = rSges(3,1) * t1216 - rSges(3,2) * t1210;
t1449 = t1250 * m(3);
t1067 = (t1161 + t1449) * t1217 - t1147 * t1211;
t1353 = t1067 * t1195;
t1418 = t1168 * t1192 + t1210 * t1267;
t1401 = pkin(3) * t1015;
t1284 = t1210 * t1401;
t1311 = t1211 * t1216;
t1324 = t1195 * t1217;
t1331 = t1195 * t1210;
t1373 = t1003 * t1224;
t1374 = t1003 * t1217;
t1276 = t1210 * t1373;
t972 = t1276 - t1401;
t919 = (((t1003 * t1324 + t1015 * t1197) * t1397 + ((-t1284 + t1373) * t1211 + pkin(2) * t1374) * t1325 + t972 * t1197) * t1003 + (t1015 * t1324 + (t1192 * t1197 - t1311 * t1331 - t1197) * t1003) * t1401) * t1055;
t1087 = pkin(3) * t1311 + t1129;
t1391 = t1003 * (-t1224 * t1284 + (t1192 * t1239 + t1216 * t1394 + t1169) * t1003);
t923 = t1055 * t1315 * t1391 + (-t1197 * t1276 + (-t1087 * t1331 + t1197 * (pkin(2) * t1216 + t1397)) * t1015) / (t1087 * t1325 + t1123) * t1015;
t935 = (-t1216 * t1391 - (pkin(2) * t1015 - t1216 * t972) * t1401) * t1055;
t1457 = -t1047 * t919 - t1075 * t923 - t1353 * t935 - 0.4e1 * t1015 * ((t1413 * t1210 + t1412 * t1216) * t1015 + ((t1210 * t1414 + t1171) * t1216 + t1411 + t1418) * t1003);
t1142 = rSges(3,1) * t1206 + rSges(3,2) * t1212;
t1069 = -t1142 * t1334 + t1252 * t1197;
t1409 = m(3) * t1069;
t998 = t1001 ^ 2;
t1456 = -t1073 * t917 - t1146 * t921 - t1409 * t933 + 0.2e1 * ((t1124 * t1206 + t1171) * t1212 + t1302 + t1420) * t998;
t1143 = rSges(3,1) * t1208 + rSges(3,2) * t1214;
t1070 = -t1143 * t1332 + t1251 * t1197;
t1408 = m(3) * t1070;
t999 = t1002 ^ 2;
t1455 = -t1074 * t918 - t1146 * t922 - t1408 * t934 + 0.2e1 * ((t1124 * t1208 + t1171) * t1214 + t1302 + t1419) * t999;
t1000 = t1003 ^ 2;
t1144 = rSges(3,1) * t1210 + rSges(3,2) * t1216;
t1071 = -t1144 * t1330 + t1250 * t1197;
t1407 = m(3) * t1071;
t1454 = -t1075 * t919 - t1146 * t923 - t1407 * t935 + 0.2e1 * t1000 * ((t1124 * t1210 + t1171) * t1216 + t1302 + t1418);
t1005 = t1009 ^ 2;
t1189 = m(1) + m(2) + m(3);
t1379 = 0.2e1 * m(3);
t1453 = -t1189 * t929 + ((-t1161 * t996 - (t996 + t1005) * t1446) * t1203 - (t1009 * t1138 * t1379 + t1147 * t997) * t1382) * t1195;
t1010 = t1013 ^ 2;
t1452 = -t1189 * t933 + ((-t1161 * t998 - (t998 + t1010) * t1447) * t1207 - (t1013 * t1142 * t1379 + t1001 * t1147) * t1378) * t1195;
t1011 = t1014 ^ 2;
t1451 = -t1189 * t934 + ((-t1161 * t999 - (t999 + t1011) * t1448) * t1209 - (t1014 * t1143 * t1379 + t1002 * t1147) * t1376) * t1195;
t1012 = t1015 ^ 2;
t1450 = -t1189 * t935 + ((-t1000 * t1161 - (t1000 + t1012) * t1449) * t1211 - (t1015 * t1144 * t1379 + t1003 * t1147) * t1374) * t1195;
t1372 = t1005 * t1138;
t1386 = t1068 * t920;
t1406 = m(3) * t1197;
t1445 = -m(3) * t1386 - t1356 * t916 - t1372 * t1406 + t1453;
t1371 = t1010 * t1142;
t1385 = t1069 * t921;
t1444 = -m(3) * t1385 - t1355 * t917 - t1371 * t1406 + t1452;
t1370 = t1011 * t1143;
t1384 = t1070 * t922;
t1443 = -m(3) * t1384 - t1354 * t918 - t1370 * t1406 + t1451;
t1369 = t1012 * t1144;
t1383 = t1071 * t923;
t1442 = -m(3) * t1383 - t1353 * t919 - t1369 * t1406 + t1450;
t1437 = t1457 * t1055;
t1436 = t1458 * t1054;
t1435 = t1459 * t1053;
t1434 = t1461 * t1052;
t1433 = t1454 * t1357;
t1432 = t1455 * t1358;
t1431 = t1456 * t1359;
t1430 = t1460 * t1360;
t1429 = t1052 * t1445;
t1428 = t1053 * t1444;
t1427 = t1054 * t1443;
t1426 = t1055 * t1442;
t1132 = pkin(2) * t1217 + t1307;
t1246 = pkin(3) * t1331 - t1129 * t1197;
t1425 = t1132 * t1196 + t1246 * t1194;
t1131 = pkin(2) * t1215 + t1308;
t1247 = pkin(3) * t1333 - t1128 * t1197;
t1424 = t1131 * t1196 + t1247 * t1194;
t1130 = pkin(2) * t1213 + t1309;
t1248 = pkin(3) * t1335 - t1127 * t1197;
t1423 = t1130 * t1196 + t1248 * t1194;
t1126 = pkin(2) * t1205 + t1310;
t1249 = pkin(3) * t1339 - t1125 * t1197;
t1422 = t1126 * t1196 + t1249 * t1194;
t1193 = t1221 ^ 2;
t1396 = m(4) * t1193;
t1344 = t1146 * t1240;
t1291 = pkin(2) * t1339;
t1290 = pkin(2) * t1335;
t1289 = pkin(2) * t1333;
t1288 = pkin(2) * t1331;
t1275 = t1072 * t1360;
t1274 = t1060 * t1344;
t1273 = t1073 * t1359;
t1272 = t1061 * t1344;
t1271 = t1074 * t1358;
t1270 = t1062 * t1344;
t1269 = t1075 * t1357;
t1268 = t1063 * t1344;
t1262 = t1360 * t1410;
t1261 = t1359 * t1409;
t1260 = t1358 * t1408;
t1259 = t1357 * t1407;
t1245 = (-t1036 * t1116 - t1037 * t1112) * t1360;
t1244 = (-t1038 * t1117 - t1041 * t1113) * t1359;
t1243 = (-t1039 * t1118 - t1042 * t1114) * t1358;
t1242 = (-t1040 * t1119 - t1043 * t1115) * t1357;
t1228 = rSges(4,1);
t1227 = rSges(4,2);
t1103 = t1194 * t1217 + t1196 * t1316;
t1102 = t1194 * t1215 + t1196 * t1318;
t1101 = t1194 * t1213 + t1196 * t1320;
t1100 = t1194 * t1316 - t1196 * t1217;
t1099 = t1194 * t1318 - t1196 * t1215;
t1098 = t1194 * t1320 - t1196 * t1213;
t1097 = t1194 * t1205 + t1196 * t1322;
t1096 = t1194 * t1322 - t1196 * t1205;
t1059 = t1132 * t1194 - t1246 * t1196;
t1058 = t1131 * t1194 - t1247 * t1196;
t1057 = t1130 * t1194 - t1248 * t1196;
t1056 = t1126 * t1194 - t1249 * t1196;
t1027 = (-t1100 * t1175 + t1103 * t1179) * t1397 + (t1059 * t1179 + t1425 * t1175) * t1216 - t1091 * t1288;
t1026 = (-t1099 * t1174 + t1102 * t1178) * t1398 + (t1058 * t1178 + t1424 * t1174) * t1214 - t1090 * t1289;
t1025 = (-t1098 * t1173 + t1101 * t1177) * t1399 + (t1057 * t1177 + t1423 * t1173) * t1212 - t1089 * t1290;
t1024 = -(t1100 * t1179 + t1103 * t1175) * t1397 + (-t1175 * t1059 + t1425 * t1179) * t1216 + t1095 * t1288;
t1023 = -(t1099 * t1178 + t1102 * t1174) * t1398 + (-t1174 * t1058 + t1424 * t1178) * t1214 + t1094 * t1289;
t1022 = -(t1098 * t1177 + t1101 * t1173) * t1399 + (-t1173 * t1057 + t1423 * t1177) * t1212 + t1093 * t1290;
t1021 = (-t1096 * t1172 + t1097 * t1176) * t1400 + (t1056 * t1176 + t1422 * t1172) * t1204 - t1088 * t1291;
t1020 = -(t1096 * t1176 + t1097 * t1172) * t1400 + (-t1172 * t1056 + t1422 * t1176) * t1204 + t1092 * t1291;
t1019 = (-t1040 * t1115 + t1043 * t1119) * t1357;
t1018 = (-t1039 * t1114 + t1042 * t1118) * t1358;
t1017 = (-t1038 * t1113 + t1041 * t1117) * t1359;
t1016 = (-t1036 * t1112 + t1037 * t1116) * t1360;
t1008 = (-t1032 * t1115 + t1035 * t1119) * t1055;
t1007 = (-t1031 * t1114 + t1034 * t1118) * t1054;
t1006 = (-t1030 * t1113 + t1033 * t1117) * t1053;
t1004 = (-t1028 * t1112 + t1029 * t1116) * t1052;
t995 = (-t1024 * t1115 + t1027 * t1119) * t1055;
t994 = (-t1023 * t1114 + t1026 * t1118) * t1054;
t993 = (-t1022 * t1113 + t1025 * t1117) * t1053;
t992 = (-t1020 * t1112 + t1021 * t1116) * t1052;
t985 = t1043 * t1259 + (t1027 * t1189 + t1035 * t1353) * t1055;
t984 = t1042 * t1260 + (t1026 * t1189 + t1034 * t1354) * t1054;
t983 = t1041 * t1261 + (t1025 * t1189 + t1033 * t1355) * t1053;
t982 = t1040 * t1259 + (t1024 * t1189 + t1032 * t1353) * t1055;
t981 = t1039 * t1260 + (t1023 * t1189 + t1031 * t1354) * t1054;
t980 = t1038 * t1261 + (t1022 * t1189 + t1030 * t1355) * t1053;
t977 = t1037 * t1262 + (t1021 * t1189 + t1029 * t1356) * t1052;
t976 = t1036 * t1262 + (t1020 * t1189 + t1028 * t1356) * t1052;
t967 = t1043 * t1269 + (t1027 * t1353 + t1035 * t1047) * t1055;
t966 = t1042 * t1271 + (t1026 * t1354 + t1034 * t1046) * t1054;
t965 = t1041 * t1273 + (t1025 * t1355 + t1033 * t1045) * t1053;
t964 = t1040 * t1269 + (t1024 * t1353 + t1032 * t1047) * t1055;
t963 = t1039 * t1271 + (t1023 * t1354 + t1031 * t1046) * t1054;
t962 = t1038 * t1273 + (t1022 * t1355 + t1030 * t1045) * t1053;
t961 = t1037 * t1275 + (t1021 * t1356 + t1029 * t1044) * t1052;
t960 = t1036 * t1275 + (t1020 * t1356 + t1028 * t1044) * t1052;
t956 = t1008 * t1353 + t1019 * t1407 + t1189 * t995;
t955 = t1007 * t1354 + t1018 * t1408 + t1189 * t994;
t954 = t1006 * t1355 + t1017 * t1409 + t1189 * t993;
t952 = t1004 * t1356 + t1016 * t1410 + t1189 * t992;
t951 = t1008 * t1047 + t1019 * t1075 + t1353 * t995;
t950 = t1007 * t1046 + t1018 * t1074 + t1354 * t994;
t949 = t1006 * t1045 + t1017 * t1073 + t1355 * t993;
t948 = t1004 * t1044 + t1016 * t1072 + t1356 * t992;
t1 = [(t1186 * t1227 - t1187 * t1228) * t1396 + t1020 * t1429 + t1022 * t1428 + t1023 * t1427 + t1024 * t1426 + t1433 * t1040 + t1432 * t1039 + t1431 * t1038 + t1430 * t1036 + t1437 * t1032 + t1436 * t1031 + t1435 * t1030 + t1434 * t1028 + ((t1036 * t1274 + (t1020 * t1410 + t1028 * t1072) * t1052) * t1245 + (-(t1020 * t976 + t1028 * t960) * t1116 - (t1021 * t976 + t1029 * t960) * t1112) * t1052 + (t1039 * t1270 + (t1023 * t1408 + t1031 * t1074) * t1054) * t1243 + (-(t1023 * t981 + t1031 * t963) * t1118 - (t1026 * t981 + t1034 * t963) * t1114) * t1054 + (t1038 * t1272 + (t1022 * t1409 + t1030 * t1073) * t1053) * t1244 + (-(t1022 * t980 + t1030 * t962) * t1117 - (t1025 * t980 + t1033 * t962) * t1113) * t1053 + (t1040 * t1268 + (t1024 * t1407 + t1032 * t1075) * t1055) * t1242 + (-(t1024 * t982 + t1032 * t964) * t1119 - (t1027 * t982 + t1035 * t964) * t1115) * t1055) * t1193; -(t1186 * t1228 + t1187 * t1227) * t1396 + t1021 * t1429 + t1025 * t1428 + t1026 * t1427 + t1027 * t1426 + t1433 * t1043 + t1432 * t1042 + t1431 * t1041 + t1430 * t1037 + t1437 * t1035 + t1436 * t1034 + t1435 * t1033 + t1434 * t1029 + ((t1037 * t1274 + (t1021 * t1410 + t1029 * t1072) * t1052) * t1245 + (-(t1020 * t977 + t1028 * t961) * t1116 - (t1021 * t977 + t1029 * t961) * t1112) * t1052 + (t1042 * t1270 + (t1026 * t1408 + t1034 * t1074) * t1054) * t1243 + (-(t1023 * t984 + t1031 * t966) * t1118 - (t1026 * t984 + t1034 * t966) * t1114) * t1054 + (t1041 * t1272 + (t1025 * t1409 + t1033 * t1073) * t1053) * t1244 + (-(t1022 * t983 + t1030 * t965) * t1117 - (t1025 * t983 + t1033 * t965) * t1113) * t1053 + (t1043 * t1268 + (t1027 * t1407 + t1035 * t1075) * t1055) * t1242 + (-(t1024 * t985 + t1032 * t967) * t1119 - (t1027 * t985 + t1035 * t967) * t1115) * t1055) * t1193; (-t1064 * t916 - t1065 * t917 - t1066 * t918 - t1067 * t919) * t1195 + (-t1386 - t1385 - t1384 - t1383 + (-t1369 - t1370 - t1371 - t1372) * t1197) * m(3) + (-t1112 * t977 - t1113 * t983 - t1114 * t984 - t1115 * t985 - t1116 * t976 - t1117 * t980 - t1118 * t981 - t1119 * t982) * t1193 + t1450 + t1451 + t1452 + t1453; t1442 * t995 + t1443 * t994 + t1444 * t993 + t1445 * t992 + t1454 * t1019 + t1455 * t1018 + t1456 * t1017 + t1460 * t1016 + t1457 * t1008 + t1458 * t1007 + t1459 * t1006 + t1461 * t1004 + ((t1004 * t1072 + t1016 * t1146 + t1410 * t992) * t1245 + (-(t1020 * t952 + t1028 * t948) * t1116 - (t1021 * t952 + t1029 * t948) * t1112) * t1052 + (t1007 * t1074 + t1018 * t1146 + t1408 * t994) * t1243 + (-(t1023 * t955 + t1031 * t950) * t1118 - (t1026 * t955 + t1034 * t950) * t1114) * t1054 + (t1006 * t1073 + t1017 * t1146 + t1409 * t993) * t1244 + (-(t1022 * t954 + t1030 * t949) * t1117 - (t1025 * t954 + t1033 * t949) * t1113) * t1053 + (t1008 * t1075 + t1019 * t1146 + t1407 * t995) * t1242 + (-(t1024 * t956 + t1032 * t951) * t1119 - (t1027 * t956 + t1035 * t951) * t1115) * t1055) * t1193;];
taucX  = t1;