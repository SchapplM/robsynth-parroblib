% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRRR8V2G3A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x12]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR8V2G3A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_regmin: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:05:08
% EndTime: 2020-08-06 18:05:21
% DurationCPUTime: 12.52s
% Computational Cost: add. (56343->455), mult. (113148->975), div. (4797->17), fcn. (103914->22), ass. (0->419)
t1104 = cos(qJ(2,3));
t1112 = pkin(7) + pkin(6);
t1068 = t1104 * t1112;
t1098 = sin(qJ(2,3));
t1046 = pkin(2) * t1098 - t1068;
t1091 = sin(pkin(4));
t1093 = cos(pkin(4));
t1097 = sin(qJ(3,3));
t1290 = t1093 * t1097;
t1025 = pkin(3) * t1290 + t1046 * t1091;
t1103 = cos(qJ(3,3));
t1304 = t1091 * t1098;
t1087 = t1103 ^ 2;
t1374 = pkin(3) * t1087;
t998 = pkin(2) * t1290 + t1025 * t1103 + t1304 * t1374;
t989 = 0.1e1 / t998;
t1106 = cos(qJ(2,2));
t1069 = t1106 * t1112;
t1100 = sin(qJ(2,2));
t1047 = pkin(2) * t1100 - t1069;
t1099 = sin(qJ(3,2));
t1288 = t1093 * t1099;
t1026 = pkin(3) * t1288 + t1047 * t1091;
t1105 = cos(qJ(3,2));
t1302 = t1091 * t1100;
t1088 = t1105 ^ 2;
t1373 = pkin(3) * t1088;
t999 = pkin(2) * t1288 + t1026 * t1105 + t1302 * t1373;
t992 = 0.1e1 / t999;
t1108 = cos(qJ(2,1));
t1070 = t1108 * t1112;
t1102 = sin(qJ(2,1));
t1048 = pkin(2) * t1102 - t1070;
t1101 = sin(qJ(3,1));
t1286 = t1093 * t1101;
t1027 = pkin(3) * t1286 + t1048 * t1091;
t1107 = cos(qJ(3,1));
t1300 = t1091 * t1102;
t1089 = t1107 ^ 2;
t1372 = pkin(3) * t1089;
t1000 = pkin(2) * t1286 + t1027 * t1107 + t1300 * t1372;
t995 = 0.1e1 / t1000;
t1394 = 0.2e1 * t1087 - 0.1e1;
t1393 = 0.2e1 * t1088 - 0.1e1;
t1392 = 0.2e1 * t1089 - 0.1e1;
t1094 = legFrame(3,2);
t1075 = sin(t1094);
t1078 = cos(t1094);
t1110 = xDP(2);
t1111 = xDP(1);
t1037 = t1075 * t1110 - t1078 * t1111;
t1090 = sin(pkin(8));
t1109 = xDP(3);
t1061 = t1109 * t1090;
t1092 = cos(pkin(8));
t1019 = t1037 * t1092 + t1061;
t1291 = t1092 * t1109;
t1022 = t1037 * t1090 - t1291;
t1289 = t1093 * t1098;
t1299 = t1091 * t1103;
t971 = t1019 * t1299 + (t1019 * t1289 + t1022 * t1104) * t1097;
t990 = 0.1e1 / t998 ^ 2;
t959 = t971 ^ 2 * t990;
t1095 = legFrame(2,2);
t1076 = sin(t1095);
t1079 = cos(t1095);
t1038 = t1076 * t1110 - t1079 * t1111;
t1020 = t1038 * t1092 + t1061;
t1023 = t1038 * t1090 - t1291;
t1287 = t1093 * t1100;
t1297 = t1091 * t1105;
t972 = t1020 * t1297 + (t1020 * t1287 + t1023 * t1106) * t1099;
t993 = 0.1e1 / t999 ^ 2;
t960 = t972 ^ 2 * t993;
t1096 = legFrame(1,2);
t1077 = sin(t1096);
t1080 = cos(t1096);
t1039 = t1077 * t1110 - t1080 * t1111;
t1021 = t1039 * t1092 + t1061;
t1024 = t1039 * t1090 - t1291;
t1285 = t1093 * t1102;
t1295 = t1091 * t1107;
t973 = t1021 * t1295 + (t1021 * t1285 + t1024 * t1108) * t1101;
t996 = 0.1e1 / t1000 ^ 2;
t961 = t973 ^ 2 * t996;
t1036 = t1090 * t1108 + t1092 * t1285;
t1018 = t1036 * t1101 + t1092 * t1295;
t1369 = pkin(3) * t1107;
t1064 = pkin(2) + t1369;
t1042 = t1064 * t1102 - t1070;
t1045 = t1064 * t1286;
t1012 = t1042 * t1295 + t1045;
t1067 = t1102 * t1112;
t979 = (t1064 * t1108 + t1067) * t1021 * t1093 - t1024 * t1042;
t1336 = 0.1e1 / t1012 * t979;
t1263 = 0.2e1 * t1336;
t1354 = t973 * t996;
t1164 = t1263 * t1354;
t1273 = t1101 * t1107;
t1360 = t995 * t961;
t1051 = pkin(2) * t1108 + t1067;
t1282 = t1093 * t1108;
t1292 = t1092 * t1093;
t985 = (t1090 * t1102 - t1092 * t1282) * t1369 - t1051 * t1292 + t1048 * t1090;
t1388 = (t1018 * t1164 + t985 * t1360) * t1273;
t1035 = t1090 * t1106 + t1092 * t1287;
t1017 = t1035 * t1099 + t1092 * t1297;
t1370 = pkin(3) * t1105;
t1063 = pkin(2) + t1370;
t1041 = t1063 * t1100 - t1069;
t1044 = t1063 * t1288;
t1011 = t1041 * t1297 + t1044;
t1066 = t1100 * t1112;
t978 = (t1063 * t1106 + t1066) * t1020 * t1093 - t1023 * t1041;
t1338 = 0.1e1 / t1011 * t978;
t1264 = 0.2e1 * t1338;
t1356 = t972 * t993;
t1165 = t1264 * t1356;
t1276 = t1099 * t1105;
t1361 = t992 * t960;
t1050 = pkin(2) * t1106 + t1066;
t1283 = t1093 * t1106;
t984 = (t1090 * t1100 - t1092 * t1283) * t1370 - t1050 * t1292 + t1047 * t1090;
t1387 = (t1017 * t1165 + t984 * t1361) * t1276;
t1034 = t1090 * t1104 + t1092 * t1289;
t1016 = t1034 * t1097 + t1092 * t1299;
t1371 = pkin(3) * t1103;
t1062 = pkin(2) + t1371;
t1040 = t1062 * t1098 - t1068;
t1043 = t1062 * t1290;
t1010 = t1040 * t1299 + t1043;
t1065 = t1098 * t1112;
t977 = (t1062 * t1104 + t1065) * t1019 * t1093 - t1022 * t1040;
t1340 = 0.1e1 / t1010 * t977;
t1265 = 0.2e1 * t1340;
t1358 = t971 * t990;
t1166 = t1265 * t1358;
t1279 = t1097 * t1103;
t1362 = t989 * t959;
t1049 = pkin(2) * t1104 + t1065;
t1284 = t1093 * t1104;
t983 = (t1090 * t1098 - t1092 * t1284) * t1371 - t1049 * t1292 + t1046 * t1090;
t1386 = (t1016 * t1166 + t983 * t1362) * t1279;
t1115 = 0.1e1 / pkin(3) ^ 2;
t1339 = 0.1e1 / t1010 ^ 2 * t977 ^ 2;
t1218 = t1115 * t1339;
t1174 = t1097 * t1218;
t1114 = 0.1e1 / pkin(3);
t1219 = t1114 * t1340;
t1175 = t1104 * t1219;
t1278 = t1098 * t1103;
t1211 = t1091 * t1278;
t1220 = t1097 * t1340;
t1359 = t971 * t989;
t1229 = t1112 * t1359;
t1318 = t1104 * t989;
t1238 = t971 * t1318;
t1256 = t989 * t1340;
t1186 = t1097 * t1229;
t953 = t1186 - t1340;
t926 = -((t1091 * t1238 + t1093 * t1219) * t1374 + ((-t1220 + t1229) * t1098 + pkin(2) * t1238) * t1299 + t1093 * t953) * t989 * t1359 - (t1091 * t1175 + (t1093 * t1087 - t1097 * t1211 - t1093) * t1359) * t1256;
t1226 = t1114 * t926 * t983;
t1028 = pkin(3) * t1278 + t1046;
t1071 = pkin(2) ^ 2 + t1112 ^ 2;
t1113 = pkin(3) ^ 2;
t1379 = 0.2e1 * pkin(2);
t1347 = pkin(3) * t1379;
t1259 = (-t1112 * t1220 + (t1087 * t1113 + t1103 * t1347 + t1071) * t1359) * t1358;
t1281 = t1093 * t1114;
t1293 = t1091 * t1114;
t1375 = pkin(2) * t1114;
t929 = -t1259 * t1281 - (-t1093 * t1186 + (-t1097 * t1028 * t1293 + t1093 * (t1103 * t1375 + t1087)) * t1340) / (t1028 * t1299 + t1043) * t1219;
t1319 = t1103 * t929;
t1385 = t1016 * (t1174 - t1319) + t1103 * t1226;
t1173 = t1103 * t1218;
t1325 = t1097 * t929;
t1384 = t1016 * (t1173 + t1325) - t1097 * t1226;
t1337 = 0.1e1 / t1011 ^ 2 * t978 ^ 2;
t1215 = t1115 * t1337;
t1171 = t1099 * t1215;
t1216 = t1114 * t1338;
t1172 = t1106 * t1216;
t1275 = t1100 * t1105;
t1210 = t1091 * t1275;
t1217 = t1099 * t1338;
t1357 = t972 * t992;
t1228 = t1112 * t1357;
t1316 = t1106 * t992;
t1234 = t972 * t1316;
t1255 = t992 * t1338;
t1183 = t1099 * t1228;
t954 = t1183 - t1338;
t927 = -((t1091 * t1234 + t1093 * t1216) * t1373 + ((-t1217 + t1228) * t1100 + pkin(2) * t1234) * t1297 + t1093 * t954) * t992 * t1357 - (t1091 * t1172 + (t1093 * t1088 - t1099 * t1210 - t1093) * t1357) * t1255;
t1225 = t1114 * t927 * t984;
t1029 = pkin(3) * t1275 + t1047;
t1258 = (-t1112 * t1217 + (t1088 * t1113 + t1105 * t1347 + t1071) * t1357) * t1356;
t930 = -t1258 * t1281 - (-t1093 * t1183 + (-t1099 * t1029 * t1293 + t1093 * (t1105 * t1375 + t1088)) * t1338) / (t1029 * t1297 + t1044) * t1216;
t1317 = t1105 * t930;
t1383 = t1017 * (t1171 - t1317) + t1105 * t1225;
t1170 = t1105 * t1215;
t1323 = t1099 * t930;
t1382 = t1017 * (t1170 + t1323) - t1099 * t1225;
t1335 = 0.1e1 / t1012 ^ 2 * t979 ^ 2;
t1212 = t1115 * t1335;
t1168 = t1101 * t1212;
t1213 = t1114 * t1336;
t1169 = t1108 * t1213;
t1272 = t1102 * t1107;
t1209 = t1091 * t1272;
t1214 = t1101 * t1336;
t1355 = t973 * t995;
t1227 = t1112 * t1355;
t1314 = t1108 * t995;
t1230 = t973 * t1314;
t1254 = t995 * t1336;
t1180 = t1101 * t1227;
t955 = t1180 - t1336;
t928 = -((t1091 * t1230 + t1093 * t1213) * t1372 + ((-t1214 + t1227) * t1102 + pkin(2) * t1230) * t1295 + t1093 * t955) * t995 * t1355 - (t1091 * t1169 + (t1093 * t1089 - t1101 * t1209 - t1093) * t1355) * t1254;
t1224 = t1114 * t928 * t985;
t1030 = pkin(3) * t1272 + t1048;
t1257 = (-t1112 * t1214 + (t1089 * t1113 + t1107 * t1347 + t1071) * t1355) * t1354;
t931 = -t1257 * t1281 - (-t1093 * t1180 + (-t1101 * t1030 * t1293 + t1093 * (t1107 * t1375 + t1089)) * t1336) / (t1030 * t1295 + t1045) * t1213;
t1315 = t1107 * t931;
t1381 = t1018 * (t1168 - t1315) + t1107 * t1224;
t1167 = t1107 * t1212;
t1321 = t1101 * t931;
t1380 = t1018 * (t1167 + t1321) - t1101 * t1224;
t1378 = pkin(2) * t1097;
t1377 = pkin(2) * t1099;
t1376 = pkin(2) * t1101;
t1368 = t926 * t989;
t1367 = t927 * t992;
t1366 = t928 * t995;
t932 = t1103 * t1259 + (pkin(2) * t1219 - t1103 * t953) * t1256;
t1365 = t932 * t989;
t933 = t1105 * t1258 + (pkin(2) * t1216 - t1105 * t954) * t1255;
t1364 = t933 * t992;
t934 = t1107 * t1257 + (pkin(2) * t1213 - t1107 * t955) * t1254;
t1363 = t934 * t995;
t1353 = t983 * t989;
t1352 = t984 * t992;
t1351 = t985 * t995;
t1306 = t1090 * t1093;
t986 = (t1090 * t1284 + t1092 * t1098) * t1371 + t1049 * t1306 + t1046 * t1092;
t1350 = t986 * t989;
t987 = (t1090 * t1283 + t1092 * t1100) * t1370 + t1050 * t1306 + t1047 * t1092;
t1349 = t987 * t992;
t988 = (t1090 * t1282 + t1092 * t1102) * t1369 + t1051 * t1306 + t1048 * t1092;
t1348 = t988 * t995;
t1142 = -0.2e1 * t1175 * t1359;
t1280 = t1097 * t1098;
t1298 = t1091 * t1104;
t920 = t1093 * t929 + t926 * t1298;
t956 = t959 + t1218;
t1346 = -t1093 * t1174 + t1103 * t920 + (t1097 * t1142 - t956 * t1278 - t929 * t1280) * t1091;
t1345 = -t1097 * t920 - t929 * t1211 + (t1103 * t1142 + t956 * t1280) * t1091 - t1093 * t1173;
t1141 = -0.2e1 * t1172 * t1357;
t1277 = t1099 * t1100;
t1296 = t1091 * t1106;
t921 = t1093 * t930 + t927 * t1296;
t957 = t960 + t1215;
t1344 = -t1093 * t1171 + t1105 * t921 + (t1099 * t1141 - t957 * t1275 - t930 * t1277) * t1091;
t1343 = -t1099 * t921 - t930 * t1210 + (t1105 * t1141 + t957 * t1277) * t1091 - t1093 * t1170;
t1140 = -0.2e1 * t1169 * t1355;
t1274 = t1101 * t1102;
t1294 = t1091 * t1108;
t922 = t1093 * t931 + t928 * t1294;
t958 = t961 + t1212;
t1342 = -t1093 * t1168 + t1107 * t922 + (t1101 * t1140 - t958 * t1272 - t931 * t1274) * t1091;
t1341 = -t1101 * t922 - t931 * t1209 + (t1107 * t1140 + t958 * t1274) * t1091 - t1093 * t1167;
t1033 = t1090 * t1285 - t1092 * t1108;
t1013 = t1033 * t1101 + t1090 * t1295;
t1334 = t1013 * t995;
t1031 = t1090 * t1289 - t1092 * t1104;
t1014 = t1031 * t1097 + t1090 * t1299;
t1333 = t1014 * t989;
t1032 = t1090 * t1287 - t1092 * t1106;
t1015 = t1032 * t1099 + t1090 * t1297;
t1332 = t1015 * t992;
t1331 = t1075 * t989;
t1330 = t1076 * t992;
t1329 = t1077 * t995;
t1328 = t1078 * t989;
t1327 = t1079 * t992;
t1326 = t1080 * t995;
t1324 = t1098 * t989;
t1322 = t1100 * t992;
t1320 = t1102 * t995;
t1313 = t1016 * t1075;
t1312 = t1016 * t1078;
t1311 = t1017 * t1076;
t1310 = t1017 * t1079;
t1309 = t1018 * t1077;
t1308 = t1018 * t1080;
t1307 = t1090 * t1091;
t1305 = t1091 * t1097;
t1303 = t1091 * t1099;
t1301 = t1091 * t1101;
t1268 = pkin(2) * t1359;
t1267 = pkin(2) * t1357;
t1266 = pkin(2) * t1355;
t1262 = t929 * t1353;
t1261 = t930 * t1352;
t1260 = t931 * t1351;
t1253 = t931 * t1334;
t1252 = t929 * t1333;
t1251 = t930 * t1332;
t1250 = t1016 * t1368;
t1249 = t1017 * t1367;
t1248 = t1018 * t1366;
t1247 = t1097 * t1362;
t1246 = t1098 * t1362;
t1245 = t1099 * t1361;
t1244 = t1100 * t1361;
t1243 = t1101 * t1360;
t1242 = t1102 * t1360;
t1241 = t1103 * t1368;
t1240 = t1103 * t1362;
t1239 = t1104 * t1362;
t1237 = t1105 * t1367;
t1236 = t1105 * t1361;
t1235 = t1106 * t1361;
t1233 = t1107 * t1366;
t1232 = t1107 * t1360;
t1231 = t1108 * t1360;
t1223 = t1097 ^ 2 * t1368;
t1222 = t1099 ^ 2 * t1367;
t1221 = t1101 ^ 2 * t1366;
t1208 = -0.2e1 * t1014 * t1340;
t1207 = t1016 * t1265;
t1206 = -0.2e1 * t1015 * t1338;
t1205 = t1017 * t1264;
t1204 = -0.2e1 * t1013 * t1336;
t1203 = t1018 * t1263;
t1202 = pkin(6) * t1219;
t1201 = pkin(6) * t1216;
t1200 = pkin(6) * t1213;
t1199 = t983 * t1247;
t1198 = t984 * t1245;
t1197 = t985 * t1243;
t1196 = t983 * t1240;
t1195 = t984 * t1236;
t1194 = t985 * t1232;
t1193 = t1333 * t1339;
t1192 = t1332 * t1337;
t1191 = t1334 * t1335;
t1190 = t1016 * t1223;
t1189 = t1017 * t1222;
t1188 = t1018 * t1221;
t1187 = t1097 * t1241;
t1184 = t1099 * t1237;
t1181 = t1101 * t1233;
t1149 = pkin(3) * t1301 - t1048 * t1093;
t982 = -t1036 * t1372 - t1051 * t1090 * t1107 + (pkin(2) * t1301 + t1149 * t1107) * t1092;
t1163 = t1013 * t934 + t928 * t982;
t1151 = pkin(3) * t1305 - t1046 * t1093;
t980 = -t1034 * t1374 - t1049 * t1090 * t1103 + (pkin(2) * t1305 + t1151 * t1103) * t1092;
t1162 = t1014 * t932 + t926 * t980;
t1150 = pkin(3) * t1303 - t1047 * t1093;
t981 = -t1035 * t1373 - t1050 * t1090 * t1105 + (pkin(2) * t1303 + t1150 * t1105) * t1092;
t1161 = t1015 * t933 + t927 * t981;
t1160 = t1016 * t1187;
t1159 = t1017 * t1184;
t1158 = t1018 * t1181;
t1157 = t1394 * t1166;
t1156 = t1393 * t1165;
t1155 = t1392 * t1164;
t1154 = t932 * t1298 + t926 * t1379;
t1153 = t933 * t1296 + t927 * t1379;
t1152 = t934 * t1294 + t928 * t1379;
t1001 = t1049 * t1092 + t1151 * t1090;
t965 = (t1031 * t1075 + t1078 * t1304) * t1374 + (-t1001 * t1075 + t1025 * t1078) * t1103 + (-t1075 * t1307 + t1078 * t1093) * t1378;
t1148 = t932 * t1313 + t926 * t965;
t962 = -(t1031 * t1078 - t1075 * t1304) * t1374 + (t1001 * t1078 + t1025 * t1075) * t1103 + (t1075 * t1093 + t1078 * t1307) * t1378;
t1147 = t932 * t1312 - t926 * t962;
t1002 = t1050 * t1092 + t1150 * t1090;
t966 = (t1032 * t1076 + t1079 * t1302) * t1373 + (-t1002 * t1076 + t1026 * t1079) * t1105 + (-t1076 * t1307 + t1079 * t1093) * t1377;
t1146 = t933 * t1311 + t927 * t966;
t963 = -(t1032 * t1079 - t1076 * t1302) * t1373 + (t1002 * t1079 + t1026 * t1076) * t1105 + (t1076 * t1093 + t1079 * t1307) * t1377;
t1145 = t933 * t1310 - t927 * t963;
t1003 = t1051 * t1092 + t1149 * t1090;
t967 = (t1033 * t1077 + t1080 * t1300) * t1372 + (-t1003 * t1077 + t1027 * t1080) * t1107 + (-t1077 * t1307 + t1080 * t1093) * t1376;
t1144 = t934 * t1309 + t928 * t967;
t964 = -(t1033 * t1080 - t1077 * t1300) * t1372 + (t1003 * t1080 + t1027 * t1077) * t1107 + (t1077 * t1093 + t1080 * t1307) * t1376;
t1143 = t934 * t1308 - t928 * t964;
t923 = pkin(6) * t926 + t932 * t1304;
t914 = t1093 * t1103 * t932 - t1097 * t923;
t947 = t1097 * t1268 + t1103 * t1202 / 0.2e1;
t1139 = t947 * t1207 + t914 * t983;
t917 = -t1103 * t923 - t932 * t1290;
t950 = t1103 * t1268 - t1097 * t1202 / 0.2e1;
t1138 = t950 * t1207 + t917 * t983;
t924 = pkin(6) * t927 + t933 * t1302;
t915 = t1093 * t1105 * t933 - t1099 * t924;
t948 = t1099 * t1267 + t1105 * t1201 / 0.2e1;
t1137 = t948 * t1205 + t915 * t984;
t918 = -t1105 * t924 - t933 * t1288;
t951 = t1105 * t1267 - t1099 * t1201 / 0.2e1;
t1136 = t951 * t1205 + t918 * t984;
t925 = pkin(6) * t928 + t934 * t1300;
t916 = t1093 * t1107 * t934 - t1101 * t925;
t949 = t1101 * t1266 + t1107 * t1200 / 0.2e1;
t1135 = t949 * t1203 + t916 * t985;
t919 = -t1107 * t925 - t934 * t1286;
t952 = t1107 * t1266 - t1101 * t1200 / 0.2e1;
t1134 = t952 * t1203 + t919 * t985;
t944 = t1394 * t959;
t1124 = t1016 * t1157 + t944 * t1353;
t945 = t1393 * t960;
t1123 = t1017 * t1156 + t945 * t1352;
t946 = t1392 * t961;
t1122 = t1018 * t1155 + t946 * t1351;
t913 = -pkin(6) * t1321 + t1152 * t1107;
t912 = -pkin(6) * t1315 - t1152 * t1101;
t911 = -pkin(6) * t1323 + t1153 * t1105;
t910 = -pkin(6) * t1317 - t1153 * t1099;
t909 = -pkin(6) * t1325 + t1154 * t1103;
t908 = -pkin(6) * t1319 - t1154 * t1097;
t1 = [t964 * t1363 + t963 * t1364 + t962 * t1365, -t1078 * t1250 - t1079 * t1249 - t1080 * t1248, (-t1143 * t1314 - t1145 * t1316 - t1147 * t1318 - t964 * t1242 - t963 * t1244 - t962 * t1246) * t1091, (t1143 * t1320 + t1145 * t1322 + t1147 * t1324 - t964 * t1231 - t963 * t1235 - t962 * t1239) * t1091, -t1078 * t1190 - t1079 * t1189 - t1080 * t1188 + (-t1078 * t1386 - t1079 * t1387 - t1080 * t1388) * t1114, -0.2e1 * t1078 * t1160 - 0.2e1 * t1079 * t1159 - 0.2e1 * t1080 * t1158 + (-t1078 * t1124 - t1079 * t1123 - t1080 * t1122) * t1114, -t1380 * t1326 - t1382 * t1327 - t1384 * t1328, t1381 * t1326 + t1383 * t1327 + t1385 * t1328, (t1078 * t1262 + t1079 * t1261 + t1080 * t1260) * t1114, (-t1308 * t913 + t1342 * t964) * t995 + (-t1310 * t911 + t1344 * t963) * t992 + (-t1312 * t909 + t1346 * t962) * t989 + (t1135 * t1326 + t1137 * t1327 + t1139 * t1328 + (t1078 * t1199 + t1079 * t1198 + t1080 * t1197) * pkin(2)) * t1114, (-t1308 * t912 + t1341 * t964) * t995 + (-t1310 * t910 + t1343 * t963) * t992 + (-t1312 * t908 + t1345 * t962) * t989 + (t1134 * t1326 + t1136 * t1327 + t1138 * t1328 + (t1078 * t1196 + t1079 * t1195 + t1080 * t1194) * pkin(2)) * t1114, 0; t967 * t1363 + t966 * t1364 + t965 * t1365, t1075 * t1250 + t1076 * t1249 + t1077 * t1248, (t1144 * t1314 + t1146 * t1316 + t1148 * t1318 - t967 * t1242 - t966 * t1244 - t965 * t1246) * t1091, (-t1144 * t1320 - t1146 * t1322 - t1148 * t1324 - t967 * t1231 - t966 * t1235 - t965 * t1239) * t1091, t1075 * t1190 + t1076 * t1189 + t1077 * t1188 + (t1075 * t1386 + t1076 * t1387 + t1077 * t1388) * t1114, 0.2e1 * t1075 * t1160 + 0.2e1 * t1076 * t1159 + 0.2e1 * t1077 * t1158 + (t1075 * t1124 + t1076 * t1123 + t1077 * t1122) * t1114, t1380 * t1329 + t1382 * t1330 + t1384 * t1331, -t1381 * t1329 - t1383 * t1330 - t1385 * t1331, (-t1075 * t1262 - t1076 * t1261 - t1077 * t1260) * t1114, (t1309 * t913 + t1342 * t967) * t995 + (t1311 * t911 + t1344 * t966) * t992 + (t1313 * t909 + t1346 * t965) * t989 + (-t1135 * t1329 - t1137 * t1330 - t1139 * t1331 + (-t1075 * t1199 - t1076 * t1198 - t1077 * t1197) * pkin(2)) * t1114, (t1309 * t912 + t1341 * t967) * t995 + (t1311 * t910 + t1343 * t966) * t992 + (t1313 * t908 + t1345 * t965) * t989 + (-t1134 * t1329 - t1136 * t1330 - t1138 * t1331 + (-t1075 * t1196 - t1076 * t1195 - t1077 * t1194) * pkin(2)) * t1114, 0; t982 * t1363 + t981 * t1364 + t980 * t1365, t927 * t1332 + t926 * t1333 + t928 * t1334, (t1161 * t1316 + t1162 * t1318 + t1163 * t1314 - t982 * t1242 - t981 * t1244 - t980 * t1246) * t1091, (-t1161 * t1322 - t1162 * t1324 - t1163 * t1320 - t982 * t1231 - t981 * t1235 - t980 * t1239) * t1091, t1013 * t1221 + t1014 * t1223 + t1015 * t1222 + ((t1013 * t1164 - t988 * t1360) * t1273 + (t1015 * t1165 - t987 * t1361) * t1276 + (t1014 * t1166 - t986 * t1362) * t1279) * t1114, 0.2e1 * t1013 * t1181 + 0.2e1 * t1014 * t1187 + 0.2e1 * t1015 * t1184 + (t1013 * t1155 + t1014 * t1157 + t1015 * t1156 - t946 * t1348 - t945 * t1349 - t944 * t1350) * t1114, t1101 * t1253 + t1097 * t1252 + t1099 * t1251 + (t1103 * t1193 + t1105 * t1192 + t1107 * t1191) * t1115 + (t1097 * t926 * t1350 + t1099 * t927 * t1349 + t1101 * t928 * t1348) * t1114, t1107 * t1253 + t1103 * t1252 + t1105 * t1251 + (-t1097 * t1193 - t1099 * t1192 - t1101 * t1191) * t1115 + (t1233 * t988 + t1237 * t987 + t1241 * t986) * t1114, (t1348 * t931 + t1349 * t930 + t1350 * t929) * t1114, (t1013 * t913 + t1342 * t982) * t995 + (t1015 * t911 + t1344 * t981) * t992 + (t1014 * t909 + t1346 * t980) * t989 + ((t1204 * t949 + t988 * t916) * t995 + (t1206 * t948 + t987 * t915) * t992 + (t1208 * t947 + t986 * t914) * t989 + (t1243 * t988 + t1245 * t987 + t1247 * t986) * pkin(2)) * t1114, (t1013 * t912 + t1341 * t982) * t995 + (t1015 * t910 + t1343 * t981) * t992 + (t1014 * t908 + t1345 * t980) * t989 + ((t1204 * t952 + t988 * t919) * t995 + (t1206 * t951 + t987 * t918) * t992 + (t1208 * t950 + t986 * t917) * t989 + (t1232 * t988 + t1236 * t987 + t1240 * t986) * pkin(2)) * t1114, 0;];
tau_reg  = t1;
