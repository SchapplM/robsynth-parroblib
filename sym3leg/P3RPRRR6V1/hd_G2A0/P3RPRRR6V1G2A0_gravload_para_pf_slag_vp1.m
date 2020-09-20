% Calculate Gravitation load for parallel robot
% P3RPRRR6V1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:37
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR6V1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:36:00
% EndTime: 2020-08-06 18:36:01
% DurationCPUTime: 1.43s
% Computational Cost: add. (762->200), mult. (936->282), div. (51->10), fcn. (765->53), ass. (0->122)
t1351 = sin(pkin(7));
t1369 = -pkin(6) - pkin(5);
t1311 = t1351 * t1369 - pkin(1);
t1357 = sin(qJ(1,3));
t1302 = t1311 * t1357;
t1352 = cos(pkin(7));
t1363 = cos(qJ(1,3));
t1403 = t1363 * t1369;
t1404 = t1363 * t1351;
t1434 = t1302 - pkin(2) * t1404 - (pkin(2) * t1357 + t1403) * t1352;
t1359 = sin(qJ(1,2));
t1303 = t1311 * t1359;
t1365 = cos(qJ(1,2));
t1401 = t1365 * t1369;
t1402 = t1365 * t1351;
t1433 = t1303 - pkin(2) * t1402 - (pkin(2) * t1359 + t1401) * t1352;
t1361 = sin(qJ(1,1));
t1304 = t1311 * t1361;
t1367 = cos(qJ(1,1));
t1399 = t1367 * t1369;
t1400 = t1367 * t1351;
t1432 = t1304 - pkin(2) * t1400 - (pkin(2) * t1361 + t1399) * t1352;
t1355 = legFrame(1,2);
t1332 = sin(t1355);
t1335 = cos(t1355);
t1301 = g(1) * t1335 - g(2) * t1332;
t1366 = cos(qJ(3,1));
t1322 = pkin(3) * t1366 + pkin(2);
t1329 = t1352 * pkin(1);
t1307 = 0.1e1 / (t1329 + t1322);
t1344 = qJ(1,1) + pkin(7);
t1325 = sin(t1344);
t1328 = cos(t1344);
t1341 = t1367 * pkin(1);
t1360 = sin(qJ(3,1));
t1378 = rSges(3,1) * t1366 - rSges(3,2) * t1360;
t1374 = pkin(2) + t1378;
t1424 = rSges(3,3) + pkin(5);
t1392 = g(3) * t1424;
t1420 = t1361 * pkin(1);
t1415 = t1307 * (-m(1) * (g(3) * (-rSges(1,1) * t1361 - rSges(1,2) * t1367) + t1301 * (rSges(1,1) * t1367 - rSges(1,2) * t1361)) - m(2) * (g(3) * (-rSges(2,1) * t1325 - rSges(2,2) * t1328 - t1420) + t1301 * (rSges(2,1) * t1328 - rSges(2,2) * t1325 + t1341)) - m(3) * (-g(3) * t1420 + t1301 * t1341 + (t1301 * t1374 + t1392) * t1328 + (-g(3) * t1374 + t1301 * t1424) * t1325));
t1431 = t1415 / 0.2e1;
t1354 = legFrame(2,2);
t1331 = sin(t1354);
t1334 = cos(t1354);
t1300 = g(1) * t1334 - g(2) * t1331;
t1364 = cos(qJ(3,2));
t1321 = pkin(3) * t1364 + pkin(2);
t1306 = 0.1e1 / (t1329 + t1321);
t1343 = qJ(1,2) + pkin(7);
t1324 = sin(t1343);
t1327 = cos(t1343);
t1339 = t1365 * pkin(1);
t1358 = sin(qJ(3,2));
t1380 = rSges(3,1) * t1364 - rSges(3,2) * t1358;
t1375 = pkin(2) + t1380;
t1421 = t1359 * pkin(1);
t1417 = t1306 * (-m(1) * (g(3) * (-rSges(1,1) * t1359 - rSges(1,2) * t1365) + t1300 * (rSges(1,1) * t1365 - rSges(1,2) * t1359)) - m(2) * (g(3) * (-rSges(2,1) * t1324 - rSges(2,2) * t1327 - t1421) + t1300 * (rSges(2,1) * t1327 - rSges(2,2) * t1324 + t1339)) - m(3) * (-g(3) * t1421 + t1300 * t1339 + (t1300 * t1375 + t1392) * t1327 + (-g(3) * t1375 + t1300 * t1424) * t1324));
t1430 = t1417 / 0.2e1;
t1353 = legFrame(3,2);
t1330 = sin(t1353);
t1333 = cos(t1353);
t1299 = g(1) * t1333 - g(2) * t1330;
t1362 = cos(qJ(3,3));
t1320 = pkin(3) * t1362 + pkin(2);
t1305 = 0.1e1 / (t1329 + t1320);
t1342 = qJ(1,3) + pkin(7);
t1323 = sin(t1342);
t1326 = cos(t1342);
t1337 = t1363 * pkin(1);
t1356 = sin(qJ(3,3));
t1382 = rSges(3,1) * t1362 - rSges(3,2) * t1356;
t1376 = pkin(2) + t1382;
t1422 = t1357 * pkin(1);
t1419 = t1305 * (-m(1) * (g(3) * (-rSges(1,1) * t1357 - rSges(1,2) * t1363) + t1299 * (rSges(1,1) * t1363 - rSges(1,2) * t1357)) - m(2) * (g(3) * (-rSges(2,1) * t1323 - rSges(2,2) * t1326 - t1422) + t1299 * (rSges(2,1) * t1326 - rSges(2,2) * t1323 + t1337)) - m(3) * (-g(3) * t1422 + t1299 * t1337 + (t1299 * t1376 + t1392) * t1326 + (-g(3) * t1376 + t1299 * t1424) * t1323));
t1429 = t1419 / 0.2e1;
t1428 = -0.2e1 * pkin(2);
t1427 = 0.2e1 * pkin(2);
t1426 = 0.2e1 * t1369;
t1425 = -m(2) - m(3);
t1423 = m(3) / pkin(3);
t1418 = t1305 / t1356;
t1416 = t1306 / t1358;
t1414 = t1307 / t1360;
t1413 = t1320 * t1357;
t1412 = t1321 * t1359;
t1411 = t1322 * t1361;
t1410 = t1356 * t1330;
t1409 = t1356 * t1333;
t1408 = t1358 * t1331;
t1407 = t1358 * t1334;
t1406 = t1360 * t1332;
t1405 = t1360 * t1335;
t1398 = -pkin(7) + qJ(3,1);
t1397 = -pkin(7) + qJ(3,2);
t1396 = -pkin(7) + qJ(3,3);
t1395 = pkin(7) + qJ(3,1);
t1394 = pkin(7) + qJ(3,2);
t1393 = pkin(7) + qJ(3,3);
t1391 = (t1352 * t1357 + t1404) * t1362 ^ 2 * pkin(3);
t1390 = (t1352 * t1359 + t1402) * t1364 ^ 2 * pkin(3);
t1389 = (t1352 * t1361 + t1400) * t1366 ^ 2 * pkin(3);
t1296 = g(1) * t1330 + g(2) * t1333;
t1388 = t1425 * t1296 * t1418;
t1297 = g(1) * t1331 + g(2) * t1334;
t1387 = t1425 * t1297 * t1416;
t1298 = g(1) * t1332 + g(2) * t1335;
t1386 = t1425 * t1298 * t1414;
t1284 = t1296 * t1382 + (g(3) * t1326 + t1299 * t1323) * (-rSges(3,1) * t1356 - rSges(3,2) * t1362);
t1385 = t1284 * ((t1403 + t1413) * t1352 - t1302 + t1320 * t1404) * t1418;
t1285 = t1297 * t1380 + (g(3) * t1327 + t1300 * t1324) * (-rSges(3,1) * t1358 - rSges(3,2) * t1364);
t1384 = t1285 * ((t1401 + t1412) * t1352 - t1303 + t1321 * t1402) * t1416;
t1286 = t1298 * t1378 + (g(3) * t1328 + t1301 * t1325) * (-rSges(3,1) * t1360 - rSges(3,2) * t1366);
t1383 = t1286 * ((t1399 + t1411) * t1352 - t1304 + t1322 * t1400) * t1414;
t1319 = t1329 + pkin(2);
t1318 = -t1355 + t1344;
t1317 = t1355 + t1344;
t1316 = -t1354 + t1343;
t1315 = t1354 + t1343;
t1314 = -t1353 + t1342;
t1313 = t1353 + t1342;
t1 = [(t1335 * t1389 + (pkin(3) * t1406 - t1335 * t1432) * t1366 + t1319 * t1406) * t1386 + (t1334 * t1390 + (pkin(3) * t1408 - t1334 * t1433) * t1364 + t1319 * t1408) * t1387 + (t1333 * t1391 + (pkin(3) * t1410 - t1333 * t1434) * t1362 + t1319 * t1410) * t1388 - m(4) * g(1) + (cos(t1318) + cos(t1317)) * t1431 + (cos(t1316) + cos(t1315)) * t1430 + (cos(t1314) + cos(t1313)) * t1429 + (t1333 * t1385 + t1334 * t1384 + t1335 * t1383) * t1423; (-t1332 * t1389 + (pkin(3) * t1405 + t1332 * t1432) * t1366 + t1319 * t1405) * t1386 + (-t1331 * t1390 + (pkin(3) * t1407 + t1331 * t1433) * t1364 + t1319 * t1407) * t1387 + (-t1330 * t1391 + (pkin(3) * t1409 + t1330 * t1434) * t1362 + t1319 * t1409) * t1388 - m(4) * g(2) + (-sin(t1317) + sin(t1318)) * t1431 + (-sin(t1315) + sin(t1316)) * t1430 + (-sin(t1313) + sin(t1314)) * t1429 + (-t1330 * t1385 - t1331 * t1384 - t1332 * t1383) * t1423; -t1325 * t1415 + ((t1322 * t1367 - t1361 * t1369) * t1352 - t1311 * t1367 - t1351 * t1411) * t1366 * t1386 - t1324 * t1417 + ((t1321 * t1365 - t1359 * t1369) * t1352 - t1311 * t1365 - t1351 * t1412) * t1364 * t1387 - t1323 * t1419 + ((t1320 * t1363 - t1357 * t1369) * t1352 - t1311 * t1363 - t1351 * t1413) * t1362 * t1388 - m(4) * g(3) + (-(t1325 * t1426 + t1328 * t1428 - 0.2e1 * t1341 + (-cos(qJ(1,1) - t1398) - cos(qJ(1,1) + t1395)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,1)) + t1360 * t1427 + (sin(t1395) + sin(t1398)) * pkin(1)) * t1286 - (t1324 * t1426 + t1327 * t1428 - 0.2e1 * t1339 + (-cos(qJ(1,2) - t1397) - cos(qJ(1,2) + t1394)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,2)) + t1358 * t1427 + (sin(t1394) + sin(t1397)) * pkin(1)) * t1285 - (t1323 * t1426 + t1326 * t1428 - 0.2e1 * t1337 + (-cos(qJ(1,3) - t1396) - cos(qJ(1,3) + t1393)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,3)) + t1356 * t1427 + (sin(t1393) + sin(t1396)) * pkin(1)) * t1284) * t1423;];
taugX  = t1;
