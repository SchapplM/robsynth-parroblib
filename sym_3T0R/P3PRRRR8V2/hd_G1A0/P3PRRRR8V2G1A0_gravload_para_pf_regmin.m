% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR8V2G1A0
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
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR8V2G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_regmin: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:36:06
% EndTime: 2020-08-06 17:36:07
% DurationCPUTime: 1.23s
% Computational Cost: add. (863->200), mult. (1825->393), div. (54->7), fcn. (2022->22), ass. (0->160)
t1337 = cos(qJ(2,1));
t1331 = sin(qJ(2,1));
t1338 = pkin(7) + pkin(6);
t1358 = t1331 * t1338;
t1301 = pkin(2) * t1337 + t1358;
t1319 = sin(pkin(8));
t1321 = cos(pkin(8));
t1309 = t1338 * t1337;
t1298 = pkin(2) * t1331 - t1309;
t1322 = cos(pkin(4));
t1320 = sin(pkin(4));
t1330 = sin(qJ(3,1));
t1361 = t1330 * t1320;
t1340 = pkin(3) * t1361 - t1298 * t1322;
t1412 = t1301 * t1321 + t1319 * t1340;
t1335 = cos(qJ(2,2));
t1329 = sin(qJ(2,2));
t1362 = t1329 * t1338;
t1300 = pkin(2) * t1335 + t1362;
t1308 = t1338 * t1335;
t1297 = pkin(2) * t1329 - t1308;
t1328 = sin(qJ(3,2));
t1365 = t1328 * t1320;
t1341 = pkin(3) * t1365 - t1297 * t1322;
t1411 = t1300 * t1321 + t1319 * t1341;
t1333 = cos(qJ(2,3));
t1327 = sin(qJ(2,3));
t1366 = t1327 * t1338;
t1299 = pkin(2) * t1333 + t1366;
t1307 = t1338 * t1333;
t1296 = pkin(2) * t1327 - t1307;
t1326 = sin(qJ(3,3));
t1369 = t1326 * t1320;
t1342 = pkin(3) * t1369 - t1296 * t1322;
t1410 = t1299 * t1321 + t1319 * t1342;
t1332 = cos(qJ(3,3));
t1409 = pkin(3) * t1332 ^ 2;
t1334 = cos(qJ(3,2));
t1408 = pkin(3) * t1334 ^ 2;
t1336 = cos(qJ(3,1));
t1407 = pkin(3) * t1336 ^ 2;
t1406 = t1320 * g(3);
t1381 = t1322 * t1326;
t1384 = t1320 * t1327;
t1257 = 0.1e1 / (t1384 * t1409 + (pkin(3) * t1381 + t1296 * t1320) * t1332 + pkin(2) * t1381);
t1302 = g(1) * t1319 - g(2) * t1321;
t1303 = g(1) * t1321 + g(2) * t1319;
t1323 = legFrame(3,3);
t1310 = sin(t1323);
t1313 = cos(t1323);
t1405 = ((-t1406 + (t1302 * t1313 + t1303 * t1310) * t1322) * t1333 + t1327 * (-t1302 * t1310 + t1303 * t1313)) * t1257;
t1379 = t1322 * t1328;
t1383 = t1320 * t1329;
t1258 = 0.1e1 / (t1383 * t1408 + (pkin(3) * t1379 + t1297 * t1320) * t1334 + pkin(2) * t1379);
t1324 = legFrame(2,3);
t1311 = sin(t1324);
t1314 = cos(t1324);
t1404 = ((-t1406 + (t1302 * t1314 + t1303 * t1311) * t1322) * t1335 + t1329 * (-t1302 * t1311 + t1303 * t1314)) * t1258;
t1377 = t1322 * t1330;
t1382 = t1320 * t1331;
t1259 = 0.1e1 / (t1382 * t1407 + (pkin(3) * t1377 + t1298 * t1320) * t1336 + pkin(2) * t1377);
t1325 = legFrame(1,3);
t1312 = sin(t1325);
t1315 = cos(t1325);
t1403 = ((-t1406 + (t1302 * t1315 + t1303 * t1312) * t1322) * t1337 + (-t1302 * t1312 + t1303 * t1315) * t1331) * t1259;
t1304 = pkin(3) * t1332 + pkin(2);
t1287 = t1304 * t1327 - t1307;
t1357 = t1332 * t1320;
t1263 = 0.1e1 / (t1287 * t1357 + t1304 * t1381);
t1266 = -t1310 * t1319 + t1313 * t1321;
t1269 = t1310 * t1321 + t1313 * t1319;
t1390 = (t1304 * t1333 + t1366) * t1322;
t1402 = (-t1266 * t1287 - t1269 * t1390) * t1263;
t1305 = pkin(3) * t1334 + pkin(2);
t1288 = t1305 * t1329 - t1308;
t1355 = t1334 * t1320;
t1264 = 0.1e1 / (t1288 * t1355 + t1305 * t1379);
t1267 = -t1311 * t1319 + t1314 * t1321;
t1270 = t1311 * t1321 + t1314 * t1319;
t1389 = (t1305 * t1335 + t1362) * t1322;
t1401 = (-t1267 * t1288 - t1270 * t1389) * t1264;
t1306 = pkin(3) * t1336 + pkin(2);
t1289 = t1306 * t1331 - t1309;
t1353 = t1336 * t1320;
t1265 = 0.1e1 / (t1289 * t1353 + t1306 * t1377);
t1268 = -t1312 * t1319 + t1315 * t1321;
t1271 = t1312 * t1321 + t1315 * t1319;
t1388 = (t1306 * t1337 + t1358) * t1322;
t1400 = (-t1268 * t1289 - t1271 * t1388) * t1265;
t1399 = (-t1266 * t1390 + t1269 * t1287) * t1263;
t1398 = (-t1267 * t1389 + t1270 * t1288) * t1264;
t1397 = (-t1268 * t1388 + t1271 * t1289) * t1265;
t1290 = g(1) * t1310 - g(2) * t1313;
t1293 = g(1) * t1313 + g(2) * t1310;
t1374 = t1322 * t1333;
t1396 = (t1293 * (t1319 * t1374 + t1321 * t1327) + t1290 * (-t1319 * t1327 + t1321 * t1374) - t1333 * t1406) * t1257;
t1380 = t1322 * t1327;
t1272 = t1319 * t1380 - t1321 * t1333;
t1275 = t1319 * t1333 + t1321 * t1380;
t1395 = (g(3) * t1384 - t1272 * t1293 - t1275 * t1290) * t1257;
t1291 = g(1) * t1311 - g(2) * t1314;
t1294 = g(1) * t1314 + g(2) * t1311;
t1372 = t1322 * t1335;
t1394 = (t1294 * (t1319 * t1372 + t1321 * t1329) + t1291 * (-t1319 * t1329 + t1321 * t1372) - t1335 * t1406) * t1258;
t1378 = t1322 * t1329;
t1273 = t1319 * t1378 - t1321 * t1335;
t1276 = t1319 * t1335 + t1321 * t1378;
t1393 = (g(3) * t1383 - t1273 * t1294 - t1276 * t1291) * t1258;
t1292 = g(1) * t1312 - g(2) * t1315;
t1295 = g(1) * t1315 + g(2) * t1312;
t1370 = t1322 * t1337;
t1392 = (t1295 * (t1319 * t1370 + t1321 * t1331) + t1292 * (-t1319 * t1331 + t1321 * t1370) - t1337 * t1406) * t1259;
t1376 = t1322 * t1331;
t1274 = t1319 * t1376 - t1321 * t1337;
t1277 = t1319 * t1337 + t1321 * t1376;
t1391 = (g(3) * t1382 - t1274 * t1295 - t1277 * t1292) * t1259;
t1375 = t1322 * t1332;
t1373 = t1322 * t1334;
t1371 = t1322 * t1336;
t1368 = t1326 * t1327;
t1367 = t1326 * t1333;
t1364 = t1328 * t1329;
t1363 = t1328 * t1335;
t1360 = t1330 * t1331;
t1359 = t1330 * t1337;
t1356 = t1332 * t1333;
t1354 = t1334 * t1335;
t1352 = t1336 * t1337;
t1351 = pkin(2) * t1369;
t1350 = pkin(2) * t1365;
t1349 = pkin(2) * t1361;
t1348 = t1326 * t1405;
t1347 = t1332 * t1405;
t1346 = t1328 * t1404;
t1345 = t1334 * t1404;
t1344 = t1330 * t1403;
t1343 = t1336 * t1403;
t1339 = 0.1e1 / pkin(3);
t1283 = t1331 * t1371 - t1361;
t1282 = t1329 * t1373 - t1365;
t1281 = t1327 * t1375 - t1369;
t1280 = t1322 * t1360 + t1353;
t1279 = t1322 * t1364 + t1355;
t1278 = t1322 * t1368 + t1357;
t1262 = t1319 * t1301 - t1321 * t1340;
t1261 = t1319 * t1300 - t1321 * t1341;
t1260 = t1319 * t1299 - t1321 * t1342;
t1244 = -t1271 * t1353 - (-t1268 * t1337 + t1271 * t1376) * t1330;
t1243 = -t1270 * t1355 - (-t1267 * t1335 + t1270 * t1378) * t1328;
t1242 = -t1269 * t1357 - (-t1266 * t1333 + t1269 * t1380) * t1326;
t1241 = -t1268 * t1353 - (t1268 * t1376 + t1271 * t1337) * t1330;
t1240 = -t1267 * t1355 - (t1267 * t1378 + t1270 * t1335) * t1328;
t1239 = -t1266 * t1357 - (t1266 * t1380 + t1269 * t1333) * t1326;
t1235 = (-t1283 * t1319 + t1321 * t1352) * t1295 - t1292 * (t1283 * t1321 + t1319 * t1352) + g(3) * (t1331 * t1353 + t1377);
t1234 = (-t1282 * t1319 + t1321 * t1354) * t1294 - t1291 * (t1282 * t1321 + t1319 * t1354) + g(3) * (t1329 * t1355 + t1379);
t1233 = (-t1281 * t1319 + t1321 * t1356) * t1293 - t1290 * (t1281 * t1321 + t1319 * t1356) + g(3) * (t1327 * t1357 + t1381);
t1232 = t1295 * (-t1280 * t1319 + t1321 * t1359) - (t1280 * t1321 + t1319 * t1359) * t1292 - g(3) * (-t1320 * t1360 + t1371);
t1231 = t1294 * (-t1279 * t1319 + t1321 * t1363) - (t1279 * t1321 + t1319 * t1363) * t1291 - g(3) * (-t1320 * t1364 + t1373);
t1230 = t1293 * (-t1278 * t1319 + t1321 * t1367) - (t1278 * t1321 + t1319 * t1367) * t1290 - g(3) * (-t1320 * t1368 + t1375);
t1 = [(-(-(t1274 * t1315 + t1277 * t1312) * t1407 + (-t1262 * t1312 + t1315 * t1412) * t1336 + t1271 * t1349) * t1259 - (-(t1273 * t1314 + t1276 * t1311) * t1408 + (-t1261 * t1311 + t1314 * t1411) * t1334 + t1270 * t1350) * t1258 - (-(t1272 * t1313 + t1275 * t1310) * t1409 + (-t1260 * t1310 + t1313 * t1410) * t1332 + t1269 * t1351) * t1257) * g(3), 0, t1239 * t1396 + t1240 * t1394 + t1241 * t1392, t1239 * t1395 + t1240 * t1393 + t1241 * t1391, 0, 0, 0, 0, 0, t1239 * t1347 + t1240 * t1345 + t1241 * t1343 + (t1230 * t1399 + t1231 * t1398 + t1232 * t1397) * t1339, -t1239 * t1348 - t1240 * t1346 - t1241 * t1344 + (t1233 * t1399 + t1234 * t1398 + t1235 * t1397) * t1339, -g(1); (-((-t1274 * t1312 + t1277 * t1315) * t1407 + (t1262 * t1315 + t1312 * t1412) * t1336 - t1268 * t1349) * t1259 - ((-t1273 * t1311 + t1276 * t1314) * t1408 + (t1261 * t1314 + t1311 * t1411) * t1334 - t1267 * t1350) * t1258 - ((-t1272 * t1310 + t1275 * t1313) * t1409 + (t1260 * t1313 + t1310 * t1410) * t1332 - t1266 * t1351) * t1257) * g(3), 0, t1242 * t1396 + t1243 * t1394 + t1244 * t1392, t1242 * t1395 + t1243 * t1393 + t1244 * t1391, 0, 0, 0, 0, 0, t1242 * t1347 + t1243 * t1345 + t1244 * t1343 + (t1230 * t1402 + t1231 * t1401 + t1232 * t1400) * t1339, -t1242 * t1348 - t1243 * t1346 - t1244 * t1344 + (t1233 * t1402 + t1234 * t1401 + t1235 * t1400) * t1339, -g(2); -0.3e1 * g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3);];
tau_reg  = t1;
