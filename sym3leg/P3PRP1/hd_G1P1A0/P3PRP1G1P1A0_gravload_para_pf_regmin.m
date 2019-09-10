% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRP1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x11]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:42
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRP1G1P1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1G1P1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1G1P1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP1G1P1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1G1P1A0_gravload_para_pf_regmin: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1G1P1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1G1P1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:42:24
% EndTime: 2019-05-03 14:42:25
% DurationCPUTime: 1.76s
% Computational Cost: add. (1394->236), mult. (2761->383), div. (72->3), fcn. (1282->14), ass. (0->181)
t1302 = (pkin(2) ^ 2);
t1372 = t1302 + 1;
t1282 = legFrame(3,3);
t1268 = sin(t1282);
t1271 = cos(t1282);
t1230 = t1268 * g(1) - t1271 * g(2);
t1233 = t1271 * g(1) + t1268 * g(2);
t1285 = sin(qJ(2,3));
t1288 = cos(qJ(2,3));
t1197 = t1285 * t1230 - t1233 * t1288;
t1283 = legFrame(2,3);
t1269 = sin(t1283);
t1272 = cos(t1283);
t1231 = t1269 * g(1) - t1272 * g(2);
t1234 = t1272 * g(1) + t1269 * g(2);
t1286 = sin(qJ(2,2));
t1289 = cos(qJ(2,2));
t1198 = t1286 * t1231 - t1234 * t1289;
t1284 = legFrame(1,3);
t1270 = sin(t1284);
t1273 = cos(t1284);
t1232 = t1270 * g(1) - t1273 * g(2);
t1235 = t1273 * g(1) + t1270 * g(2);
t1287 = sin(qJ(2,1));
t1290 = cos(qJ(2,1));
t1199 = t1287 * t1232 - t1235 * t1290;
t1357 = t1290 * qJ(3,1);
t1328 = -0.2e1 * t1357;
t1371 = qJ(3,1) * t1270;
t1214 = t1273 * t1328 + t1287 * (pkin(2) * t1273 + t1371);
t1295 = qJ(3,1) ^ 2;
t1264 = -t1295 + t1372;
t1281 = t1290 ^ 2;
t1327 = 0.2e1 * t1357;
t1220 = 0.1e1 / (t1287 * pkin(2) * t1327 + t1264 * t1281 - t1295 - t1372);
t1348 = t1214 * t1220;
t1358 = t1289 * qJ(3,2);
t1330 = -0.2e1 * t1358;
t1367 = qJ(3,2) * t1269;
t1212 = t1272 * t1330 + t1286 * (pkin(2) * t1272 + t1367);
t1294 = qJ(3,2) ^ 2;
t1263 = -t1294 + t1372;
t1280 = t1289 ^ 2;
t1329 = 0.2e1 * t1358;
t1219 = 0.1e1 / (t1286 * pkin(2) * t1329 + t1263 * t1280 - t1294 - t1372);
t1350 = t1212 * t1219;
t1359 = t1288 * qJ(3,3);
t1332 = -0.2e1 * t1359;
t1363 = qJ(3,3) * t1268;
t1210 = t1271 * t1332 + t1285 * (pkin(2) * t1271 + t1363);
t1293 = qJ(3,3) ^ 2;
t1262 = -t1293 + t1372;
t1279 = t1288 ^ 2;
t1331 = 0.2e1 * t1359;
t1218 = 0.1e1 / (t1285 * pkin(2) * t1331 + t1262 * t1279 - t1293 - t1372);
t1352 = t1210 * t1218;
t1379 = t1197 * t1352 + t1198 * t1350 + t1199 * t1348;
t1370 = qJ(3,1) * t1273;
t1213 = t1270 * t1328 + t1287 * (pkin(2) * t1270 - t1370);
t1349 = t1213 * t1220;
t1366 = qJ(3,2) * t1272;
t1211 = t1269 * t1330 + t1286 * (pkin(2) * t1269 - t1366);
t1351 = t1211 * t1219;
t1362 = qJ(3,3) * t1271;
t1209 = t1268 * t1332 + t1285 * (pkin(2) * t1268 - t1362);
t1353 = t1209 * t1218;
t1378 = t1197 * t1353 + t1198 * t1351 + t1199 * t1349;
t1292 = xP(3);
t1274 = sin(t1292);
t1275 = cos(t1292);
t1298 = koppelP(1,2);
t1301 = koppelP(1,1);
t1369 = qJ(3,1) * t1298;
t1254 = pkin(2) * t1301 - t1369;
t1368 = qJ(3,1) * t1301;
t1255 = pkin(2) * t1298 + t1368;
t1374 = (t1274 * t1254 + t1255 * t1275) * t1273 - (t1254 * t1275 - t1274 * t1255) * t1270;
t1181 = -t1374 * t1287 + ((t1274 * t1301 + t1275 * t1298) * t1273 - t1270 * (-t1274 * t1298 + t1275 * t1301)) * t1327;
t1354 = t1181 * t1220;
t1297 = koppelP(2,2);
t1300 = koppelP(2,1);
t1365 = qJ(3,2) * t1297;
t1252 = pkin(2) * t1300 - t1365;
t1364 = qJ(3,2) * t1300;
t1253 = pkin(2) * t1297 + t1364;
t1375 = (t1274 * t1252 + t1253 * t1275) * t1272 - (t1252 * t1275 - t1274 * t1253) * t1269;
t1180 = -t1375 * t1286 + ((t1274 * t1300 + t1275 * t1297) * t1272 - t1269 * (-t1274 * t1297 + t1275 * t1300)) * t1329;
t1355 = t1180 * t1219;
t1296 = koppelP(3,2);
t1299 = koppelP(3,1);
t1361 = qJ(3,3) * t1296;
t1250 = pkin(2) * t1299 - t1361;
t1360 = qJ(3,3) * t1299;
t1251 = pkin(2) * t1296 + t1360;
t1376 = (t1274 * t1250 + t1251 * t1275) * t1271 - (t1250 * t1275 - t1274 * t1251) * t1268;
t1179 = -t1376 * t1285 + ((t1274 * t1299 + t1275 * t1296) * t1271 - t1268 * (-t1274 * t1296 + t1275 * t1299)) * t1331;
t1356 = t1179 * t1218;
t1377 = t1197 * t1356 + t1198 * t1355 + t1199 * t1354;
t1373 = pkin(2) * g(2);
t1344 = t1218 * t1230;
t1343 = t1219 * t1231;
t1342 = t1220 * t1232;
t1341 = t1285 * t1288;
t1340 = t1286 * t1289;
t1339 = t1287 * t1290;
t1338 = t1293 * t1296;
t1337 = t1294 * t1297;
t1336 = t1295 * t1298;
t1335 = t1372 * t1296;
t1334 = t1372 * t1297;
t1333 = t1372 * t1298;
t1326 = pkin(2) * t1368;
t1325 = pkin(2) * t1364;
t1324 = pkin(2) * t1360;
t1323 = pkin(2) * t1371;
t1322 = pkin(2) * t1367;
t1321 = pkin(2) * t1363;
t1320 = pkin(2) * t1362;
t1319 = pkin(2) * t1366;
t1318 = pkin(2) * t1370;
t1265 = pkin(2) * t1361;
t1224 = t1262 * t1299 - 0.2e1 * t1265;
t1225 = 0.2e1 * t1324 + t1335 - t1338;
t1191 = t1224 * t1275 - t1274 * t1225;
t1192 = t1274 * t1224 + t1225 * t1275;
t1266 = pkin(2) * t1365;
t1226 = t1263 * t1300 - 0.2e1 * t1266;
t1227 = 0.2e1 * t1325 + t1334 - t1337;
t1193 = t1226 * t1275 - t1274 * t1227;
t1194 = t1274 * t1226 + t1227 * t1275;
t1267 = pkin(2) * t1369;
t1228 = t1264 * t1301 - 0.2e1 * t1267;
t1229 = 0.2e1 * t1326 + t1333 - t1336;
t1195 = t1228 * t1275 - t1274 * t1229;
t1196 = t1274 * t1228 + t1229 * t1275;
t1244 = t1302 * t1299 - t1265 + t1299;
t1245 = t1324 + t1335;
t1246 = t1302 * t1300 - t1266 + t1300;
t1247 = t1325 + t1334;
t1248 = t1302 * t1301 - t1267 + t1301;
t1249 = t1326 + t1333;
t1317 = ((t1195 * t1273 + t1196 * t1270) * t1281 + (-t1195 * t1270 + t1196 * t1273) * t1339 + (-t1248 * t1275 + t1274 * t1249) * t1273 - (t1274 * t1248 + t1249 * t1275) * t1270) * t1342 + ((t1193 * t1272 + t1194 * t1269) * t1280 + (-t1193 * t1269 + t1194 * t1272) * t1340 + (-t1246 * t1275 + t1274 * t1247) * t1272 - (t1274 * t1246 + t1247 * t1275) * t1269) * t1343 + ((t1191 * t1271 + t1192 * t1268) * t1279 + (-t1191 * t1268 + t1192 * t1271) * t1341 + (-t1244 * t1275 + t1274 * t1245) * t1271 - (t1274 * t1244 + t1245 * t1275) * t1268) * t1344;
t1221 = t1262 * t1271 + 0.2e1 * t1321;
t1222 = t1263 * t1272 + 0.2e1 * t1322;
t1223 = t1264 * t1273 + 0.2e1 * t1323;
t1305 = -t1268 * t1262 + 0.2e1 * t1320;
t1308 = -t1269 * t1263 + 0.2e1 * t1319;
t1311 = -t1270 * t1264 + 0.2e1 * t1318;
t1316 = (t1223 * t1281 - t1273 * t1302 + t1311 * t1339 - t1273 - t1323) * t1342 + (t1222 * t1280 - t1272 * t1302 + t1308 * t1340 - t1272 - t1322) * t1343 + (t1221 * t1279 - t1271 * t1302 + t1305 * t1341 - t1271 - t1321) * t1344;
t1315 = (-t1223 * t1339 + t1302 * t1270 + t1311 * t1281 + t1270 - t1318) * t1342 + (-t1222 * t1340 + t1302 * t1269 + t1308 * t1280 + t1269 - t1319) * t1343 + (-t1221 * t1341 + t1302 * t1268 + t1305 * t1279 + t1268 - t1320) * t1344;
t1200 = t1230 * t1288 + t1233 * t1285;
t1201 = t1231 * t1289 + t1234 * t1286;
t1202 = t1232 * t1290 + t1235 * t1287;
t1310 = -t1295 * t1270 - t1318;
t1309 = t1273 * t1295 - t1323;
t1307 = -t1294 * t1269 - t1319;
t1306 = t1272 * t1294 - t1322;
t1304 = -t1293 * t1268 - t1320;
t1303 = t1271 * t1293 - t1321;
t1291 = pkin(2) * g(1);
t1261 = g(1) * qJ(3,1) + t1373;
t1260 = -g(2) * qJ(3,1) + t1291;
t1259 = g(1) * qJ(3,2) + t1373;
t1258 = -g(2) * qJ(3,2) + t1291;
t1257 = g(1) * qJ(3,3) + t1373;
t1256 = -g(2) * qJ(3,3) + t1291;
t1243 = t1295 * t1301 + t1267 + t1301;
t1242 = -t1298 + t1326 - t1336;
t1241 = t1294 * t1300 + t1266 + t1300;
t1240 = -t1297 + t1325 - t1337;
t1239 = t1293 * t1299 + t1265 + t1299;
t1238 = -t1296 + t1324 - t1338;
t1237 = t1275 * g(1) + t1274 * g(2);
t1236 = t1274 * g(1) - t1275 * g(2);
t1190 = (t1287 * t1260 - t1261 * t1290) * t1273 + (t1260 * t1290 + t1287 * t1261) * t1270;
t1189 = (t1286 * t1258 - t1259 * t1289) * t1272 + (t1258 * t1289 + t1286 * t1259) * t1269;
t1188 = (t1285 * t1256 - t1257 * t1288) * t1271 + (t1256 * t1288 + t1285 * t1257) * t1268;
t1178 = t1200 * t1352 + t1201 * t1350 + t1202 * t1348;
t1177 = t1200 * t1353 + t1201 * t1351 + t1202 * t1349;
t1173 = t1200 * t1356 + t1201 * t1355 + t1202 * t1354;
t1 = [t1315, 0, t1178, -t1379, t1178, t1379, (t1214 * t1190 - (t1310 * t1290 + t1287 * (-t1273 - t1309)) * t1202) * t1220 + (t1212 * t1189 - (t1307 * t1289 + t1286 * (-t1272 - t1306)) * t1201) * t1219 + (t1210 * t1188 - (t1304 * t1288 + t1285 * (-t1271 - t1303)) * t1200) * t1218 + t1315, 0, 0, 0, -t1274 * t1236 - t1275 * t1237; t1316, 0, t1177, -t1378, t1177, t1378, (t1213 * t1190 - (t1309 * t1290 - t1287 * (t1270 - t1310)) * t1202) * t1220 + (t1211 * t1189 - (t1306 * t1289 - t1286 * (t1269 - t1307)) * t1201) * t1219 + (t1209 * t1188 - (t1303 * t1288 - t1285 * (t1268 - t1304)) * t1200) * t1218 + t1316, 0, 0, 0, t1275 * t1236 - t1274 * t1237; t1317, 0, t1173, -t1377, t1173, t1377, (t1181 * t1190 - (((-t1242 * t1275 + t1274 * t1243) * t1273 - (t1274 * t1242 + t1243 * t1275) * t1270) * t1287 + t1374 * t1357) * t1202) * t1220 + (t1180 * t1189 - (((-t1240 * t1275 + t1274 * t1241) * t1272 - (t1274 * t1240 + t1241 * t1275) * t1269) * t1286 + t1375 * t1358) * t1201) * t1219 + (t1179 * t1188 - (((-t1238 * t1275 + t1274 * t1239) * t1271 - (t1274 * t1238 + t1239 * t1275) * t1268) * t1285 + t1376 * t1359) * t1200) * t1218 + t1317, 0, t1236, t1237, 0;];
tau_reg  = t1;
