% Calculate Gravitation load for parallel robot
% P3PRRRR8V2G4A0
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
% koppelP [3x3]
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
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V2G4A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G4A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:14:33
% EndTime: 2020-08-06 18:14:36
% DurationCPUTime: 2.83s
% Computational Cost: add. (1941->263), mult. (4404->515), div. (36->10), fcn. (4494->34), ass. (0->197)
t1345 = sin(qJ(2,3));
t1351 = cos(qJ(2,3));
t1335 = legFrame(3,3);
t1306 = sin(t1335);
t1338 = legFrame(3,1);
t1309 = sin(t1338);
t1312 = cos(t1335);
t1315 = cos(t1338);
t1341 = legFrame(3,2);
t1318 = sin(t1341);
t1400 = t1318 * t1315;
t1403 = t1309 * t1318;
t1321 = cos(t1341);
t1429 = g(1) * t1321;
t1230 = -t1306 * t1429 + (-t1306 * t1403 + t1312 * t1315) * g(2) + (t1306 * t1400 + t1309 * t1312) * g(3);
t1231 = t1312 * t1429 + (t1306 * t1315 + t1312 * t1403) * g(2) + (t1306 * t1309 - t1312 * t1400) * g(3);
t1331 = sin(pkin(8));
t1333 = cos(pkin(8));
t1376 = t1230 * t1331 + t1231 * t1333;
t1263 = g(1) * t1318 + (-g(2) * t1309 + g(3) * t1315) * t1321;
t1332 = sin(pkin(4));
t1334 = cos(pkin(4));
t1377 = t1230 * t1333 - t1231 * t1331;
t1438 = t1263 * t1332 + t1334 * t1377;
t1362 = t1438 * t1345 + t1376 * t1351;
t1347 = sin(qJ(2,2));
t1353 = cos(qJ(2,2));
t1336 = legFrame(2,3);
t1307 = sin(t1336);
t1339 = legFrame(2,1);
t1310 = sin(t1339);
t1313 = cos(t1336);
t1316 = cos(t1339);
t1342 = legFrame(2,2);
t1319 = sin(t1342);
t1399 = t1319 * t1316;
t1402 = t1310 * t1319;
t1322 = cos(t1342);
t1428 = g(1) * t1322;
t1232 = -t1307 * t1428 + (-t1307 * t1402 + t1313 * t1316) * g(2) + (t1307 * t1399 + t1310 * t1313) * g(3);
t1233 = t1313 * t1428 + (t1307 * t1316 + t1313 * t1402) * g(2) + (t1307 * t1310 - t1313 * t1399) * g(3);
t1374 = t1232 * t1331 + t1233 * t1333;
t1264 = g(1) * t1319 + (-g(2) * t1310 + g(3) * t1316) * t1322;
t1375 = t1232 * t1333 - t1233 * t1331;
t1439 = t1264 * t1332 + t1334 * t1375;
t1361 = t1439 * t1347 + t1374 * t1353;
t1349 = sin(qJ(2,1));
t1355 = cos(qJ(2,1));
t1337 = legFrame(1,3);
t1308 = sin(t1337);
t1340 = legFrame(1,1);
t1311 = sin(t1340);
t1314 = cos(t1337);
t1317 = cos(t1340);
t1343 = legFrame(1,2);
t1320 = sin(t1343);
t1398 = t1320 * t1317;
t1401 = t1311 * t1320;
t1323 = cos(t1343);
t1427 = g(1) * t1323;
t1234 = -t1308 * t1427 + (-t1308 * t1401 + t1314 * t1317) * g(2) + (t1308 * t1398 + t1311 * t1314) * g(3);
t1235 = t1314 * t1427 + (t1308 * t1317 + t1314 * t1401) * g(2) + (t1308 * t1311 - t1314 * t1398) * g(3);
t1372 = t1234 * t1331 + t1235 * t1333;
t1265 = g(1) * t1320 + (-g(2) * t1311 + g(3) * t1317) * t1323;
t1373 = t1234 * t1333 - t1235 * t1331;
t1440 = t1265 * t1332 + t1334 * t1373;
t1360 = t1440 * t1349 + t1372 * t1355;
t1272 = t1306 * t1333 + t1312 * t1331;
t1446 = t1272 * t1332;
t1273 = t1307 * t1333 + t1313 * t1331;
t1445 = t1273 * t1332;
t1274 = t1308 * t1333 + t1314 * t1331;
t1444 = t1274 * t1332;
t1443 = -t1265 * t1334 + t1332 * t1373;
t1442 = -t1264 * t1334 + t1332 * t1375;
t1441 = -t1263 * t1334 + t1332 * t1377;
t1437 = m(3) / pkin(3);
t1350 = cos(qJ(3,3));
t1436 = pkin(3) * t1350 ^ 2;
t1352 = cos(qJ(3,2));
t1435 = pkin(3) * t1352 ^ 2;
t1354 = cos(qJ(3,1));
t1434 = pkin(3) * t1354 ^ 2;
t1433 = pkin(3) * t1332;
t1432 = pkin(3) * t1350;
t1431 = pkin(3) * t1352;
t1430 = pkin(3) * t1354;
t1344 = sin(qJ(3,3));
t1426 = t1344 * pkin(2);
t1346 = sin(qJ(3,2));
t1425 = t1346 * pkin(2);
t1348 = sin(qJ(3,1));
t1424 = t1348 * pkin(2);
t1300 = pkin(2) + t1432;
t1357 = pkin(7) + pkin(6);
t1303 = t1357 * t1351;
t1287 = t1300 * t1345 - t1303;
t1393 = t1334 * t1344;
t1290 = t1300 * t1393;
t1396 = t1332 * t1350;
t1423 = ((t1441 * rSges(3,1) + t1362 * rSges(3,2)) * t1350 + t1344 * (t1362 * rSges(3,1) - t1441 * rSges(3,2))) / (t1287 * t1396 + t1290);
t1301 = pkin(2) + t1431;
t1304 = t1357 * t1353;
t1288 = t1301 * t1347 - t1304;
t1391 = t1334 * t1346;
t1291 = t1301 * t1391;
t1395 = t1332 * t1352;
t1422 = ((t1442 * rSges(3,1) + t1361 * rSges(3,2)) * t1352 + t1346 * (t1361 * rSges(3,1) - t1442 * rSges(3,2))) / (t1288 * t1395 + t1291);
t1302 = pkin(2) + t1430;
t1305 = t1357 * t1355;
t1289 = t1302 * t1349 - t1305;
t1389 = t1334 * t1348;
t1292 = t1302 * t1389;
t1394 = t1332 * t1354;
t1421 = ((t1443 * rSges(3,1) + t1360 * rSges(3,2)) * t1354 + t1348 * (t1360 * rSges(3,1) - t1443 * rSges(3,2))) / (t1289 * t1394 + t1292);
t1299 = (-pkin(6) - rSges(3,3)) * m(3) + m(2) * rSges(2,2);
t1381 = m(2) * rSges(2,1) + pkin(2) * m(3);
t1215 = t1362 * t1299 + (t1345 * t1376 - t1351 * t1438) * ((rSges(3,1) * t1350 - rSges(3,2) * t1344) * m(3) + t1381);
t1293 = pkin(2) * t1345 - t1303;
t1420 = t1215 / (t1290 + (t1345 * t1432 + t1293) * t1396);
t1216 = t1361 * t1299 + (t1347 * t1374 - t1353 * t1439) * ((rSges(3,1) * t1352 - rSges(3,2) * t1346) * m(3) + t1381);
t1294 = pkin(2) * t1347 - t1304;
t1419 = t1216 / (t1291 + (t1347 * t1431 + t1294) * t1395);
t1217 = t1360 * t1299 + (t1349 * t1372 - t1355 * t1440) * ((rSges(3,1) * t1354 - rSges(3,2) * t1348) * m(3) + t1381);
t1295 = pkin(2) * t1349 - t1305;
t1418 = t1217 / (t1292 + (t1349 * t1430 + t1295) * t1394);
t1266 = pkin(3) * t1393 + t1293 * t1332;
t1387 = t1345 * t1332;
t1242 = 0.1e1 / (pkin(2) * t1393 + t1266 * t1350 + t1387 * t1436);
t1417 = t1242 * t1263;
t1267 = pkin(3) * t1391 + t1294 * t1332;
t1385 = t1347 * t1332;
t1243 = 0.1e1 / (pkin(2) * t1391 + t1267 * t1352 + t1385 * t1435);
t1416 = t1243 * t1264;
t1268 = pkin(3) * t1389 + t1295 * t1332;
t1383 = t1349 * t1332;
t1244 = 0.1e1 / (pkin(2) * t1389 + t1268 * t1354 + t1383 * t1434);
t1415 = t1244 * t1265;
t1269 = -t1306 * t1331 + t1312 * t1333;
t1408 = t1269 * t1332;
t1270 = -t1307 * t1331 + t1313 * t1333;
t1407 = t1270 * t1332;
t1386 = t1345 * t1357;
t1406 = (t1300 * t1351 + t1386) * t1334;
t1384 = t1347 * t1357;
t1405 = (t1301 * t1353 + t1384) * t1334;
t1382 = t1349 * t1357;
t1404 = (t1302 * t1355 + t1382) * t1334;
t1271 = -t1308 * t1331 + t1314 * t1333;
t1397 = t1332 * t1271;
t1392 = t1334 * t1345;
t1390 = t1334 * t1347;
t1388 = t1334 * t1349;
t1278 = t1331 * t1392 - t1333 * t1351;
t1281 = t1331 * t1351 + t1333 * t1392;
t1371 = t1278 * t1312 + t1281 * t1306;
t1279 = t1331 * t1390 - t1333 * t1353;
t1282 = t1331 * t1353 + t1333 * t1390;
t1370 = t1279 * t1313 + t1282 * t1307;
t1280 = t1331 * t1388 - t1333 * t1355;
t1283 = t1331 * t1355 + t1333 * t1388;
t1369 = t1280 * t1314 + t1283 * t1308;
t1368 = -t1293 * t1334 + t1344 * t1433;
t1367 = -t1294 * t1334 + t1346 * t1433;
t1366 = -t1295 * t1334 + t1348 * t1433;
t1327 = m(1) + m(2) + m(3);
t1298 = pkin(2) * t1355 + t1382;
t1297 = pkin(2) * t1353 + t1384;
t1296 = pkin(2) * t1351 + t1386;
t1256 = -t1320 * t1444 + t1323 * t1334;
t1255 = -t1319 * t1445 + t1322 * t1334;
t1254 = -t1318 * t1446 + t1321 * t1334;
t1253 = t1298 * t1331 - t1333 * t1366;
t1252 = t1297 * t1331 - t1333 * t1367;
t1251 = t1296 * t1331 - t1333 * t1368;
t1250 = -t1298 * t1333 - t1331 * t1366;
t1249 = -t1297 * t1333 - t1331 * t1367;
t1248 = -t1296 * t1333 - t1331 * t1368;
t1247 = -t1280 * t1308 + t1283 * t1314;
t1246 = -t1279 * t1307 + t1282 * t1313;
t1245 = -t1278 * t1306 + t1281 * t1312;
t1241 = t1271 * t1401 + t1274 * t1317;
t1240 = t1270 * t1402 + t1273 * t1316;
t1239 = t1269 * t1403 + t1272 * t1315;
t1238 = t1271 * t1398 - t1274 * t1311;
t1237 = t1270 * t1399 - t1273 * t1310;
t1236 = t1269 * t1400 - t1272 * t1309;
t1229 = t1320 * t1369 + t1323 * t1383;
t1228 = t1319 * t1370 + t1322 * t1385;
t1227 = t1318 * t1371 + t1321 * t1387;
t1226 = -t1250 * t1308 + t1253 * t1314;
t1225 = -t1249 * t1307 + t1252 * t1313;
t1224 = -t1248 * t1306 + t1251 * t1312;
t1220 = t1268 * t1323 + (t1250 * t1314 + t1253 * t1308) * t1320;
t1219 = t1267 * t1322 + (t1249 * t1313 + t1252 * t1307) * t1319;
t1218 = t1266 * t1321 + (t1248 * t1312 + t1251 * t1306) * t1318;
t1 = [-t1323 * (t1271 * t1394 + t1348 * (t1271 * t1388 + t1274 * t1355)) * t1244 * t1217 - t1322 * (t1270 * t1395 + t1346 * (t1270 * t1390 + t1273 * t1353)) * t1243 * t1216 - t1321 * (t1269 * t1396 + t1344 * (t1269 * t1392 + t1272 * t1351)) * t1242 * t1215 - m(4) * g(1) + (-(-((-t1271 * t1355 + t1274 * t1388) * t1323 - t1320 * t1383) * t1434 + ((t1271 * t1298 + t1274 * t1366) * t1323 + t1268 * t1320) * t1354 + (t1320 * t1334 + t1323 * t1444) * t1424) * t1415 - (-((-t1270 * t1353 + t1273 * t1390) * t1322 - t1319 * t1385) * t1435 + ((t1270 * t1297 + t1273 * t1367) * t1322 + t1267 * t1319) * t1352 + (t1319 * t1334 + t1322 * t1445) * t1425) * t1416 - (-((-t1269 * t1351 + t1272 * t1392) * t1321 - t1318 * t1387) * t1436 + ((t1269 * t1296 + t1272 * t1368) * t1321 + t1266 * t1318) * t1350 + (t1318 * t1334 + t1321 * t1446) * t1426) * t1417) * t1327 + (t1323 * (-t1271 * t1404 + t1274 * t1289) * t1421 + t1322 * (-t1270 * t1405 + t1273 * t1288) * t1422 + t1321 * (-t1269 * t1406 + t1272 * t1287) * t1423) * t1437; ((-t1247 * t1401 - t1317 * t1369) * t1348 - t1241 * t1394) * t1418 + ((-t1246 * t1402 - t1316 * t1370) * t1346 - t1240 * t1395) * t1419 + ((-t1245 * t1403 - t1315 * t1371) * t1344 - t1239 * t1396) * t1420 - m(4) * g(2) + (-(-(t1229 * t1311 - t1247 * t1317) * t1434 + (-t1220 * t1311 + t1226 * t1317) * t1354 - (t1256 * t1311 + t1317 * t1397) * t1424) * t1415 - (-(t1228 * t1310 - t1246 * t1316) * t1435 + (-t1219 * t1310 + t1225 * t1316) * t1352 - (t1255 * t1310 + t1316 * t1407) * t1425) * t1416 - (-(t1227 * t1309 - t1245 * t1315) * t1436 + (-t1218 * t1309 + t1224 * t1315) * t1350 - (t1254 * t1309 + t1315 * t1408) * t1426) * t1417) * t1327 + ((-t1241 * t1404 + (-t1271 * t1317 + t1274 * t1401) * t1289) * t1421 + (-t1240 * t1405 + (-t1270 * t1316 + t1273 * t1402) * t1288) * t1422 + (-t1239 * t1406 + (-t1269 * t1315 + t1272 * t1403) * t1287) * t1423) * t1437; ((t1247 * t1398 - t1311 * t1369) * t1348 + t1238 * t1394) * t1418 + ((t1246 * t1399 - t1310 * t1370) * t1346 + t1237 * t1395) * t1419 + ((t1245 * t1400 - t1309 * t1371) * t1344 + t1236 * t1396) * t1420 - m(4) * g(3) + (-((t1229 * t1317 + t1247 * t1311) * t1434 + (t1220 * t1317 + t1226 * t1311) * t1354 + (t1256 * t1317 - t1311 * t1397) * t1424) * t1415 - ((t1228 * t1316 + t1246 * t1310) * t1435 + (t1219 * t1316 + t1225 * t1310) * t1352 + (t1255 * t1316 - t1310 * t1407) * t1425) * t1416 - ((t1227 * t1315 + t1245 * t1309) * t1436 + (t1218 * t1315 + t1224 * t1309) * t1350 + (t1254 * t1315 - t1309 * t1408) * t1426) * t1417) * t1327 + ((t1238 * t1404 - (t1271 * t1311 + t1274 * t1398) * t1289) * t1421 + (t1237 * t1405 - (t1270 * t1310 + t1273 * t1399) * t1288) * t1422 + (t1236 * t1406 - (t1269 * t1309 + t1272 * t1400) * t1287) * t1423) * t1437;];
taugX  = t1;
