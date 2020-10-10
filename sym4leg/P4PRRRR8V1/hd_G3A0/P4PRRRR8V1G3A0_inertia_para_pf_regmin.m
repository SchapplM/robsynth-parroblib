% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P4PRRRR8V1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [4*4x15]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-09-20 23:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4PRRRR8V1G3A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G3A0_inertia_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G3A0_inertia_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G3A0_inertia_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G3A0_inertia_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G3A0_inertia_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-09-20 23:00:17
% EndTime: 2020-09-20 23:01:13
% DurationCPUTime: 59.48s
% Computational Cost: add. (19406->932), mult. (60770->2078), div. (5040->18), fcn. (64370->30), ass. (0->846)
t868 = cos(qJ(3,4));
t851 = 0.1e1 / t868;
t880 = cos(qJ(3,3));
t856 = 0.1e1 / t880;
t882 = cos(qJ(3,2));
t858 = 0.1e1 / t882;
t884 = cos(qJ(3,1));
t860 = 0.1e1 / t884;
t895 = 0.1e1 / pkin(2);
t865 = cos(pkin(3));
t878 = sin(qJ(3,1));
t1424 = t865 * t878;
t879 = sin(qJ(2,1));
t1407 = t879 * t884;
t885 = cos(qJ(2,1));
t823 = pkin(2) * t1407 - pkin(5) * t885;
t863 = sin(pkin(3));
t1295 = pkin(2) * t1424 + t823 * t863;
t1558 = 0.1e1 / t1295 ^ 2;
t873 = legFrame(1,2);
t847 = cos(t873);
t1448 = t1558 * t847;
t1420 = t865 * t885;
t1423 = t865 * t879;
t1522 = pkin(2) * t884;
t862 = sin(pkin(6));
t864 = cos(pkin(6));
t753 = pkin(5) * (t864 * t1423 + t862 * t885) + (t864 * t1420 - t862 * t879) * t1522;
t756 = (t862 * t1420 + t864 * t879) * t1522 + (t862 * t1423 - t864 * t885) * pkin(5);
t1408 = t878 * t885;
t1409 = t878 * t879;
t948 = t865 * t1409 + t863 * t884;
t773 = -t864 * t1408 + t862 * t948;
t776 = t862 * t1408 + t864 * t948;
t1568 = t753 * t773 + t756 * t776;
t905 = t1568 * t1448;
t876 = sin(qJ(3,2));
t1426 = t865 * t876;
t877 = sin(qJ(2,2));
t1410 = t877 * t882;
t883 = cos(qJ(2,2));
t822 = pkin(2) * t1410 - pkin(5) * t883;
t1296 = pkin(2) * t1426 + t822 * t863;
t1556 = 0.1e1 / t1296 ^ 2;
t872 = legFrame(2,2);
t846 = cos(t872);
t1454 = t1556 * t846;
t1421 = t865 * t883;
t1425 = t865 * t877;
t1523 = pkin(2) * t882;
t752 = pkin(5) * (t864 * t1425 + t862 * t883) + (t864 * t1421 - t862 * t877) * t1523;
t755 = (t862 * t1421 + t864 * t877) * t1523 + (t862 * t1425 - t864 * t883) * pkin(5);
t1411 = t876 * t883;
t1412 = t876 * t877;
t949 = t865 * t1412 + t863 * t882;
t772 = -t864 * t1411 + t862 * t949;
t775 = t862 * t1411 + t864 * t949;
t1567 = t752 * t772 + t755 * t775;
t908 = t1567 * t1454;
t874 = sin(qJ(3,3));
t1428 = t865 * t874;
t875 = sin(qJ(2,3));
t1413 = t875 * t880;
t881 = cos(qJ(2,3));
t821 = pkin(2) * t1413 - pkin(5) * t881;
t1297 = pkin(2) * t1428 + t821 * t863;
t1554 = 0.1e1 / t1297 ^ 2;
t871 = legFrame(3,2);
t845 = cos(t871);
t1460 = t1554 * t845;
t1422 = t865 * t881;
t1427 = t865 * t875;
t1524 = pkin(2) * t880;
t751 = pkin(5) * (t864 * t1427 + t862 * t881) + (t864 * t1422 - t862 * t875) * t1524;
t754 = (t862 * t1422 + t864 * t875) * t1524 + (t862 * t1427 - t864 * t881) * pkin(5);
t1414 = t874 * t881;
t1415 = t874 * t875;
t950 = t865 * t1415 + t863 * t880;
t771 = -t864 * t1414 + t862 * t950;
t774 = t862 * t1414 + t864 * t950;
t1566 = t751 * t771 + t754 * t774;
t911 = t1566 * t1460;
t869 = cos(qJ(2,4));
t1429 = t865 * t869;
t867 = sin(qJ(2,4));
t1430 = t865 * t867;
t1525 = pkin(2) * t868;
t746 = -pkin(5) * (t864 * t1430 + t862 * t869) + (-t864 * t1429 + t862 * t867) * t1525;
t866 = sin(qJ(3,4));
t1417 = t866 * t869;
t1418 = t866 * t867;
t951 = t865 * t1418 + t863 * t868;
t763 = -t864 * t1417 + t862 * t951;
t1493 = t746 * t763;
t1431 = t865 * t866;
t1416 = t867 * t868;
t819 = pkin(2) * t1416 - pkin(5) * t869;
t1298 = pkin(2) * t1431 + t819 * t863;
t1552 = 0.1e1 / t1298 ^ 2;
t747 = (t862 * t1429 + t864 * t867) * t1525 + (t862 * t1430 - t864 * t869) * pkin(5);
t764 = t862 * t1417 + t864 * t951;
t1578 = t1552 * (t747 * t764 - t1493);
t870 = legFrame(4,2);
t844 = cos(t870);
t915 = t844 * t1578;
t1588 = (-t851 * t915 - t856 * t911 - t858 * t908 - t860 * t905) * t895;
t861 = 0.1e1 / t884 ^ 2;
t1536 = t861 * t878;
t859 = 0.1e1 / t882 ^ 2;
t1537 = t859 * t876;
t857 = 0.1e1 / t880 ^ 2;
t1538 = t857 * t874;
t852 = 0.1e1 / t868 ^ 2;
t1539 = t852 * t866;
t1587 = (-t1536 * t905 - t1537 * t908 - t1538 * t911 - t1539 * t915) * t895;
t840 = sin(t870);
t1444 = t840 * t895;
t1465 = t1552 * t851;
t1193 = t1444 * t1465;
t1086 = t764 * t1193;
t1377 = t746 * t764 * t1552;
t1260 = t844 * t1377;
t1445 = t840 * t851;
t1492 = t746 * t844;
t1350 = t776 * t1448;
t1402 = -0.2e1 * t753;
t843 = sin(t873);
t1579 = t1402 * t843 * t1350;
t1354 = t774 * t1460;
t1401 = -0.2e1 * t751;
t841 = sin(t871);
t1580 = t1401 * t841 * t1354;
t1352 = t775 * t1454;
t1400 = -0.2e1 * t752;
t842 = sin(t872);
t1581 = t1400 * t842 * t1352;
t1586 = t1086 * t1492 + (t1260 * t1445 + t1579 * t860 + t1580 * t856 + t1581 * t858) * t895;
t1464 = t1552 * t852;
t1344 = t844 * t1464;
t1194 = t840 * t1344;
t1126 = t746 * t1194;
t1341 = t866 * t1464;
t1055 = t1341 * t1444;
t979 = t764 * t1055;
t1585 = t979 * t1492 + (t764 * t866 * t1126 + t1536 * t1579 + t1537 * t1581 + t1538 * t1580) * t895;
t1584 = t751 * t754;
t1583 = t752 * t755;
t1582 = t753 * t756;
t1446 = t1558 * t861;
t1311 = t847 * t1446;
t1176 = t843 * t1311;
t1452 = t1556 * t859;
t1322 = t846 * t1452;
t1182 = t842 * t1322;
t1458 = t1554 * t857;
t1333 = t845 * t1458;
t1188 = t841 * t1333;
t896 = 0.1e1 / pkin(2) ^ 2;
t1572 = t746 * t896;
t1574 = t753 ^ 2;
t1575 = t752 ^ 2;
t1576 = t751 ^ 2;
t1577 = -t1126 * t1572 + (-t1176 * t1574 - t1182 * t1575 - t1188 * t1576) * t896;
t1573 = t746 * t747;
t1571 = (-t1311 * t1582 - t1322 * t1583 - t1333 * t1584 + t1344 * t1573) * t896;
t1570 = t746 ^ 2;
t1551 = 0.1e1 / t1298;
t886 = xP(4);
t848 = sin(t886);
t849 = cos(t886);
t887 = koppelP(4,2);
t891 = koppelP(4,1);
t811 = t848 * t891 + t849 * t887;
t815 = -t848 * t887 + t849 * t891;
t1047 = t811 * t844 + t815 * t840;
t1469 = t1551 * t851;
t916 = t1047 * t1469;
t685 = t764 * t916;
t1521 = t685 * t1551;
t1555 = 0.1e1 / t1296;
t889 = koppelP(2,2);
t893 = koppelP(2,1);
t813 = t848 * t893 + t849 * t889;
t817 = -t848 * t889 + t849 * t893;
t1045 = t813 * t846 + t817 * t842;
t1457 = t1555 * t858;
t910 = t1045 * t1457;
t686 = t775 * t910;
t1520 = t686 * t1555;
t1557 = 0.1e1 / t1295;
t890 = koppelP(1,2);
t894 = koppelP(1,1);
t814 = t848 * t894 + t849 * t890;
t818 = -t848 * t890 + t849 * t894;
t1044 = t814 * t847 + t818 * t843;
t1451 = t1557 * t860;
t907 = t1044 * t1451;
t687 = t776 * t907;
t1519 = t687 * t1557;
t1553 = 0.1e1 / t1297;
t888 = koppelP(3,2);
t892 = koppelP(3,1);
t812 = t848 * t892 + t849 * t888;
t816 = -t848 * t888 + t849 * t892;
t1046 = t812 * t845 + t816 * t841;
t1463 = t1553 * t856;
t913 = t1046 * t1463;
t688 = t774 * t913;
t1518 = t688 * t1553;
t1496 = t1047 * t1551;
t1308 = t878 * t1446;
t1319 = t876 * t1452;
t1330 = t874 * t1458;
t1559 = t841 * t895 * t1566;
t1560 = t842 * t895 * t1567;
t1561 = t843 * t895 * t1568;
t1564 = -t1055 * t1493 + t1308 * t1561 + t1319 * t1560 + t1330 * t1559 + t747 * t979;
t1447 = t1558 * t860;
t1453 = t1556 * t858;
t1459 = t1554 * t856;
t1563 = t747 * t1086 - t1193 * t1493 + t1447 * t1561 + t1453 * t1560 + t1459 * t1559;
t1245 = t1446 * t1582;
t1248 = t1452 * t1583;
t1251 = t1458 * t1584;
t1562 = -t1464 * t747 * t840 * t1572 + (t1245 * t843 + t1248 * t842 + t1251 * t841) * t896;
t1547 = t747 * t763;
t1546 = t754 * t771;
t1545 = t755 * t772;
t1544 = t756 * t773;
t1543 = t763 * t764;
t1542 = t771 * t774;
t1541 = t772 * t775;
t1540 = t773 * t776;
t1478 = t774 * t1553;
t1380 = t1046 * t1478;
t1271 = t856 * t1380;
t1526 = pkin(2) * t863;
t1031 = t874 * t1526 - t821 * t865;
t824 = pkin(5) * t875 + t881 * t1524;
t731 = t1031 * t862 + t864 * t824;
t723 = t1297 * t841 + t731 * t845;
t724 = t1297 * t845 - t731 * t841;
t656 = (-t723 * t812 + t724 * t816) * t1553;
t1535 = (t688 + t1271) * t656;
t1472 = t776 * t1557;
t1382 = t1044 * t1472;
t1277 = t860 * t1382;
t1029 = t878 * t1526 - t823 * t865;
t826 = pkin(5) * t879 + t885 * t1522;
t733 = t1029 * t862 + t864 * t826;
t727 = t1295 * t843 + t733 * t847;
t728 = t1295 * t847 - t733 * t843;
t658 = (-t727 * t814 + t728 * t818) * t1557;
t1534 = (t687 + t1277) * t658;
t1475 = t775 * t1555;
t1379 = t1045 * t1475;
t1268 = t858 * t1379;
t1030 = t876 * t1526 - t822 * t865;
t825 = pkin(5) * t877 + t883 * t1523;
t732 = t1030 * t862 + t864 * t825;
t725 = t1296 * t842 + t732 * t846;
t726 = t1296 * t846 - t732 * t842;
t657 = (-t725 * t813 + t726 * t817) * t1555;
t1533 = (t686 + t1268) * t657;
t1481 = t764 * t1551;
t1381 = t1047 * t1481;
t1274 = t851 * t1381;
t1032 = t866 * t1526 - t819 * t865;
t820 = pkin(5) * t867 + t869 * t1525;
t729 = t1032 * t862 + t864 * t820;
t721 = t1298 * t840 + t729 * t844;
t722 = t1298 * t844 - t729 * t840;
t653 = (-t721 * t811 + t722 * t815) * t1551;
t1532 = (t685 + t1274) * t653;
t1531 = t1044 * t1446;
t1530 = t1045 * t1452;
t1529 = t1046 * t1458;
t1528 = 0.2e1 * t863;
t1527 = 0.2e1 * t895;
t1517 = t721 * t1551;
t1516 = t722 * t1551;
t1515 = t723 * t1553;
t1514 = t724 * t1553;
t1513 = t725 * t1555;
t1512 = t726 * t1555;
t1511 = t727 * t1557;
t1510 = t728 * t1557;
t730 = t1032 * t864 - t862 * t820;
t1509 = t730 * t764;
t1508 = t730 * t1551;
t1507 = t730 * t1552;
t734 = t1031 * t864 - t862 * t824;
t1506 = t734 * t774;
t1505 = t734 * t1553;
t1504 = t734 * t1554;
t735 = t1030 * t864 - t862 * t825;
t1503 = t735 * t775;
t1502 = t735 * t1555;
t1501 = t735 * t1556;
t736 = t1029 * t864 - t862 * t826;
t1500 = t736 * t776;
t1499 = t736 * t1557;
t1498 = t736 * t1558;
t1497 = t1044 * t1557;
t1495 = t1046 * t1553;
t1494 = t1045 * t1555;
t1489 = t751 * t1553;
t1487 = t752 * t1555;
t1485 = t753 * t1557;
t1479 = t764 * t840;
t1477 = t774 * t1554;
t1476 = t774 * t841;
t1474 = t775 * t1556;
t1473 = t775 * t842;
t1471 = t776 * t1558;
t1470 = t776 * t843;
t1468 = t1551 * t867;
t1467 = t1551 * t869;
t1462 = t1553 * t875;
t1461 = t1553 * t881;
t1456 = t1555 * t877;
t1455 = t1555 * t883;
t1450 = t1557 * t879;
t1449 = t1557 * t885;
t1440 = t851 * t866;
t1439 = t856 * t874;
t1438 = t858 * t876;
t1437 = t860 * t878;
t1436 = t863 * t869;
t1435 = t863 * t881;
t1434 = t863 * t883;
t1433 = t863 * t885;
t1432 = t863 * t895;
t1419 = t865 * t895;
t1406 = 0.2e1 * t1044;
t1405 = 0.2e1 * t1047;
t1404 = 0.2e1 * t1046;
t1403 = 0.2e1 * t1045;
t681 = t895 * t746 * t916;
t1399 = t681 * t1481;
t682 = t895 * t751 * t913;
t1398 = t682 * t1478;
t683 = t895 * t752 * t910;
t1397 = t683 * t1475;
t684 = t895 * t753 * t907;
t1396 = t684 * t1472;
t1395 = t685 * t1517;
t1394 = t746 * t1521;
t1393 = t866 * t1521;
t1392 = t869 * t1521;
t1391 = t686 * t1513;
t1390 = t876 * t1520;
t1389 = t883 * t1520;
t1388 = t687 * t1511;
t1387 = t878 * t1519;
t1386 = t885 * t1519;
t1385 = t688 * t1515;
t1384 = t874 * t1518;
t1383 = t881 * t1518;
t1378 = t746 * t1496;
t1376 = t747 * t1469;
t1375 = t1046 * t1489;
t1374 = t751 * t1477;
t1373 = t1045 * t1487;
t1372 = t752 * t1474;
t1371 = t1044 * t1485;
t1370 = t753 * t1471;
t1369 = t754 * t1463;
t1368 = t755 * t1457;
t1367 = t756 * t1451;
t1366 = t1047 * t1464;
t762 = t764 ^ 2;
t1365 = t762 * t1464;
t1364 = t763 * t1469;
t1363 = t763 * t1467;
t1362 = t764 * t1467;
t1361 = t771 * t1463;
t1360 = t771 * t1461;
t1359 = t772 * t1457;
t1358 = t772 * t1455;
t1357 = t773 * t1451;
t1356 = t773 * t1449;
t1355 = t774 * t1461;
t1353 = t775 * t1455;
t1351 = t776 * t1449;
t1349 = t1551 * t1445;
t1348 = t844 * t1469;
t1347 = t1551 * t1440;
t1346 = t851 * t1468;
t1345 = t851 * t1467;
t850 = t866 ^ 2;
t1343 = t850 * t1464;
t1342 = t1552 * t1440;
t1340 = t841 * t1463;
t1339 = t845 * t1463;
t1338 = t1553 * t1439;
t1337 = t856 * t1462;
t1336 = t856 * t1461;
t833 = t841 ^ 2;
t1335 = t833 * t1458;
t837 = t845 ^ 2;
t1334 = t837 * t1458;
t853 = t874 ^ 2;
t1332 = t853 * t1458;
t1331 = t1554 * t1439;
t1329 = t842 * t1457;
t1328 = t846 * t1457;
t1327 = t1555 * t1438;
t1326 = t858 * t1456;
t1325 = t858 * t1455;
t834 = t842 ^ 2;
t1324 = t834 * t1452;
t838 = t846 ^ 2;
t1323 = t838 * t1452;
t854 = t876 ^ 2;
t1321 = t854 * t1452;
t1320 = t1556 * t1438;
t1318 = t843 * t1451;
t1317 = t847 * t1451;
t1316 = t1557 * t1437;
t1315 = t860 * t1450;
t1314 = t860 * t1449;
t835 = t843 ^ 2;
t1313 = t835 * t1446;
t839 = t847 ^ 2;
t1312 = t839 * t1446;
t855 = t878 ^ 2;
t1310 = t855 * t1446;
t1309 = t1558 * t1437;
t1307 = t863 * t1418;
t1306 = t863 * t1416;
t1305 = t863 * t1415;
t1304 = t863 * t1413;
t1303 = t863 * t1412;
t1302 = t863 * t1410;
t1301 = t863 * t1409;
t1300 = t863 * t1407;
t1294 = t681 * t1347;
t1293 = t682 * t1338;
t1292 = t683 * t1327;
t1291 = t684 * t1316;
t1290 = t764 * t1393;
t1289 = t685 * t850 * t1469;
t1288 = t685 * t1347;
t1287 = t775 * t1390;
t1286 = t686 * t854 * t1457;
t1285 = t686 * t1327;
t1284 = t776 * t1387;
t1283 = t687 * t855 * t1451;
t1282 = t687 * t1316;
t1281 = t774 * t1384;
t1280 = t688 * t853 * t1463;
t1279 = t688 * t1338;
t770 = t776 ^ 2;
t1278 = t770 * t1531;
t1276 = t1044 * t1351;
t1275 = t1047 * t1365;
t1273 = t1047 * t1362;
t768 = t774 ^ 2;
t1272 = t768 * t1529;
t1270 = t1046 * t1355;
t769 = t775 ^ 2;
t1269 = t769 * t1530;
t1267 = t1045 * t1353;
t1266 = t753 * t1317;
t1265 = t751 * t1339;
t1264 = t752 * t1328;
t1263 = t1570 * t1366;
t1262 = t1570 * t1464;
t1261 = t851 * t1378;
t1259 = t746 * t1349;
t1258 = t746 * t1348;
t803 = t865 * t868 - t1307;
t1257 = t803 * t1376;
t804 = t1306 + t1431;
t1256 = t804 * t1376;
t1254 = t1576 * t1529;
t1253 = t1575 * t1530;
t1252 = t1574 * t1531;
t1250 = t856 * t1375;
t1249 = t751 * t1340;
t1247 = t858 * t1373;
t1246 = t752 * t1329;
t1244 = t860 * t1371;
t1243 = t753 * t1318;
t805 = t865 * t880 - t1305;
t1242 = t805 * t1369;
t806 = t1304 + t1428;
t1241 = t806 * t1369;
t807 = t865 * t882 - t1303;
t1239 = t807 * t1368;
t808 = t1302 + t1426;
t1238 = t808 * t1368;
t809 = t865 * t884 - t1301;
t1236 = t809 * t1367;
t810 = t1300 + t1424;
t1235 = t810 * t1367;
t1233 = t762 * t1343;
t1232 = t762 * t1342;
t1231 = t1464 * t1543;
t1230 = t1551 * t1363;
t1229 = t764 * t1349;
t1228 = t840 * t1362;
t1227 = t764 * t1348;
t1226 = t844 * t1362;
t1225 = t768 * t1332;
t1224 = t768 * t1331;
t1223 = t769 * t1321;
t1222 = t769 * t1320;
t1221 = t770 * t1310;
t1220 = t770 * t1309;
t1219 = t1458 * t1542;
t1218 = t1553 * t1360;
t1217 = t1452 * t1541;
t1216 = t1555 * t1358;
t1215 = t1446 * t1540;
t1214 = t1557 * t1356;
t1213 = t774 * t1340;
t1212 = t841 * t1355;
t1211 = t774 * t1339;
t1210 = t845 * t1355;
t1208 = t775 * t1329;
t1207 = t842 * t1353;
t1206 = t775 * t1328;
t1205 = t846 * t1353;
t1203 = t776 * t1318;
t1202 = t843 * t1351;
t1201 = t776 * t1317;
t1200 = t847 * t1351;
t1198 = t1551 * t1346;
t1197 = t1551 * t1345;
t1196 = t866 * t1345;
t1195 = t1432 * t1468;
t1192 = t1553 * t1337;
t1191 = t1553 * t1336;
t1190 = t874 * t1336;
t1189 = t1432 * t1462;
t1186 = t1555 * t1326;
t1185 = t1555 * t1325;
t1184 = t876 * t1325;
t1183 = t1432 * t1456;
t1180 = t1557 * t1315;
t1179 = t1557 * t1314;
t1178 = t878 * t1314;
t1177 = t1432 * t1450;
t1160 = t658 * t863 * t1314;
t1162 = t657 * t863 * t1325;
t1164 = t656 * t863 * t1336;
t1166 = t653 * t863 * t1345;
t1174 = t1160 * t1470 + t1162 * t1473 + t1164 * t1476 + t1166 * t1479;
t1161 = t658 * t1201;
t1163 = t657 * t1206;
t1165 = t656 * t1211;
t1167 = t653 * t1227;
t1173 = (t1161 * t879 + t1163 * t877 + t1165 * t875 + t1167 * t867) * t863;
t1172 = t773 * t1160 + t772 * t1162 + t771 * t1164 + t763 * t1166;
t1171 = -0.2e1 * t1377;
t1170 = 0.2e1 * t833 * t1374;
t1169 = 0.2e1 * t834 * t1372;
t1168 = 0.2e1 * t835 * t1370;
t1159 = t764 * t1294;
t1158 = t774 * t1293;
t1157 = t775 * t1292;
t1156 = t776 * t1291;
t1155 = t746 * t1288;
t1154 = t764 * t1289;
t1153 = t752 * t1285;
t1152 = t775 * t1286;
t1151 = t753 * t1282;
t1150 = t776 * t1283;
t1149 = t751 * t1279;
t1148 = t774 * t1280;
t1147 = t730 * t1230;
t1146 = t734 * t1218;
t1145 = t735 * t1216;
t1144 = t736 * t1214;
t1143 = t843 * t1278;
t1142 = t1557 * t1276;
t1141 = t840 * t1275;
t1140 = t1551 * t1273;
t1139 = t841 * t1272;
t1138 = t1553 * t1270;
t1137 = t842 * t1269;
t1136 = t1555 * t1267;
t1135 = t809 * t1266;
t1134 = t810 * t1266;
t1133 = t805 * t1265;
t1132 = t806 * t1265;
t1131 = t807 * t1264;
t1130 = t808 * t1264;
t1129 = t803 * t1261;
t1128 = t804 * t1261;
t1127 = t746 * t1195;
t1125 = t803 * t1259;
t1124 = t803 * t1258;
t1123 = t804 * t1259;
t1122 = t804 * t1258;
t1121 = t1551 * t1257;
t1120 = t1551 * t1256;
t1118 = t805 * t1250;
t1117 = t806 * t1250;
t1116 = t805 * t1249;
t1115 = t806 * t1249;
t1114 = t751 * t1189;
t1113 = t807 * t1247;
t1112 = t808 * t1247;
t1111 = t807 * t1246;
t1110 = t808 * t1246;
t1109 = t752 * t1183;
t1108 = t809 * t1244;
t1107 = t810 * t1244;
t1106 = t809 * t1243;
t1105 = t810 * t1243;
t1104 = t753 * t1177;
t1103 = t1553 * t1242;
t1102 = t1553 * t1241;
t1101 = t1555 * t1239;
t1100 = t1555 * t1238;
t1099 = t1557 * t1236;
t1098 = t1557 * t1235;
t1097 = t844 * t1233;
t1096 = t844 * t1232;
t1095 = t850 * t1231;
t1094 = t1342 * t1543;
t1093 = t763 * t1196;
t1091 = t1551 * t1228;
t1090 = t1551 * t1226;
t1089 = t764 * t1198;
t1088 = t764 * t1197;
t1087 = t764 * t1196;
t1085 = t845 * t1225;
t1084 = t845 * t1224;
t1083 = t846 * t1223;
t1082 = t846 * t1222;
t1081 = t847 * t1221;
t1080 = t847 * t1220;
t1079 = t853 * t1219;
t1078 = t1331 * t1542;
t1077 = t771 * t1190;
t1076 = t854 * t1217;
t1075 = t1320 * t1541;
t1074 = t772 * t1184;
t1073 = t855 * t1215;
t1072 = t1309 * t1540;
t1071 = t773 * t1178;
t1070 = t1553 * t1212;
t1069 = t1553 * t1210;
t1068 = t774 * t1192;
t1067 = t774 * t1191;
t1066 = t774 * t1190;
t1065 = t1555 * t1207;
t1064 = t1555 * t1205;
t1063 = t775 * t1186;
t1062 = t775 * t1185;
t1061 = t775 * t1184;
t1060 = t1557 * t1202;
t1059 = t1557 * t1200;
t1058 = t776 * t1180;
t1057 = t776 * t1179;
t1056 = t776 * t1178;
t1051 = t721 * t840 - t722 * t844;
t1050 = t723 * t841 - t724 * t845;
t1049 = t725 * t842 - t726 * t846;
t1048 = t727 * t843 - t728 * t847;
t1043 = t851 * t1171;
t1042 = t837 * t1401 * t1477;
t1041 = t838 * t1400 * t1474;
t1040 = t839 * t1402 * t1471;
t1028 = t721 * t1090;
t1027 = t721 * t1089;
t1026 = t722 * t1091;
t1025 = t723 * t1069;
t1024 = t723 * t1068;
t1023 = t724 * t1070;
t1022 = t725 * t1064;
t1021 = t725 * t1063;
t1020 = t726 * t1065;
t1019 = t727 * t1059;
t1018 = t727 * t1058;
t1017 = t728 * t1060;
t1016 = t1044 * t1057;
t1015 = t1044 * t1056;
t1014 = t1047 * t1088;
t1013 = t1047 * t1087;
t1012 = t1046 * t1067;
t1011 = t1046 * t1066;
t1010 = t1045 * t1062;
t1009 = t1045 * t1061;
t1008 = t1557 * t1135;
t1007 = t1557 * t1134;
t1006 = t1553 * t1133;
t1005 = t1553 * t1132;
t1004 = t1555 * t1131;
t1003 = t1555 * t1130;
t1002 = t1551 * t1129;
t1001 = t1551 * t1128;
t999 = t1551 * t1125;
t998 = t1551 * t1124;
t997 = t1551 * t1123;
t996 = t1551 * t1122;
t995 = t1553 * t1118;
t994 = t1553 * t1117;
t993 = t1553 * t1116;
t992 = t1553 * t1115;
t991 = t1555 * t1113;
t990 = t1555 * t1112;
t989 = t1555 * t1111;
t988 = t1555 * t1110;
t987 = t1557 * t1108;
t986 = t1557 * t1107;
t985 = t1557 * t1106;
t984 = t1557 * t1105;
t983 = t1551 * t1093;
t981 = t840 * t1087;
t980 = t844 * t1087;
t978 = t1553 * t1077;
t977 = t1555 * t1074;
t976 = t1557 * t1071;
t975 = t841 * t1066;
t974 = t845 * t1066;
t973 = t842 * t1061;
t972 = t846 * t1061;
t971 = t843 * t1056;
t970 = t847 * t1056;
t969 = t1195 * t1440;
t968 = t1189 * t1439;
t967 = t1183 * t1438;
t966 = t1177 * t1437;
t965 = t844 * t1509 - t721 * t763;
t964 = t730 * t1479 + t722 * t763;
t963 = t845 * t1506 - t723 * t771;
t962 = t734 * t1476 + t724 * t771;
t961 = t846 * t1503 - t725 * t772;
t960 = t735 * t1473 + t726 * t772;
t959 = t847 * t1500 - t727 * t773;
t958 = t736 * t1470 + t728 * t773;
t947 = t840 * t1405 * t1377;
t946 = t1405 * t1260;
t945 = t841 * t1404 * t1374;
t944 = t751 * t1404 * t1354;
t943 = t842 * t1403 * t1372;
t942 = t752 * t1403 * t1352;
t941 = t843 * t1406 * t1370;
t940 = t753 * t1406 * t1350;
t939 = t1171 * t1539;
t935 = t1557 * t1015;
t934 = t1551 * t1013;
t933 = t1553 * t1011;
t932 = t1555 * t1009;
t931 = t746 * t969;
t930 = t751 * t968;
t929 = t752 * t967;
t928 = t753 * t966;
t927 = t1551 * t981;
t926 = t1551 * t980;
t925 = t1553 * t975;
t924 = t1553 * t974;
t923 = t1555 * t973;
t922 = t1555 * t972;
t921 = t1557 * t971;
t920 = t1557 * t970;
t904 = t1047 * t1578;
t903 = t1554 * t1566 * t1046;
t902 = t1556 * t1567 * t1045;
t901 = t1558 * t1568 * t1044;
t900 = (-t1419 * t746 + t764 * t1436) * t1469;
t899 = (t751 * t1419 + t774 * t1435) * t1463;
t898 = (t752 * t1419 + t775 * t1434) * t1457;
t897 = (t753 * t1419 + t776 * t1433) * t1451;
t836 = t844 ^ 2;
t832 = t840 ^ 2;
t827 = t848 ^ 2 + t849 ^ 2;
t767 = t773 ^ 2;
t766 = t772 ^ 2;
t765 = t771 ^ 2;
t761 = t763 ^ 2;
t700 = (t756 * t1419 + t773 * t1433) * t1451;
t699 = (t755 * t1419 + t772 * t1434) * t1457;
t698 = (t754 * t1419 + t771 * t1435) * t1463;
t697 = (t747 * t1419 + t763 * t1436) * t1469;
t696 = t847 * t897;
t695 = t843 * t897;
t694 = t846 * t898;
t693 = t842 * t898;
t692 = t845 * t899;
t691 = t841 * t899;
t690 = t844 * t900;
t689 = t840 * t900;
t680 = -t1177 * t756 - t700 * t878;
t679 = -t1183 * t755 - t699 * t876;
t678 = -t1189 * t754 - t698 * t874;
t677 = t700 * t884 - t756 * t966;
t676 = t699 * t882 - t755 * t967;
t675 = t698 * t880 - t754 * t968;
t674 = -t1195 * t747 - t697 * t866;
t673 = t697 * t868 - t747 * t969;
t672 = t1104 * t847 + t696 * t878;
t671 = -t1104 * t843 - t695 * t878;
t670 = t1109 * t846 + t694 * t876;
t669 = -t1109 * t842 - t693 * t876;
t668 = t1114 * t845 + t692 * t874;
t667 = -t1114 * t841 - t691 * t874;
t666 = -t696 * t884 + t847 * t928;
t665 = t695 * t884 - t843 * t928;
t664 = -t694 * t882 + t846 * t929;
t663 = t693 * t882 - t842 * t929;
t662 = -t692 * t880 + t845 * t930;
t661 = t691 * t880 - t841 * t930;
t660 = -t1127 * t844 + t690 * t866;
t659 = t1127 * t840 - t689 * t866;
t655 = -t690 * t868 - t844 * t931;
t654 = t689 * t868 + t840 * t931;
t652 = -t1176 * t770 - t1182 * t769 - t1188 * t768 - t1194 * t762;
t651 = -t1081 * t843 - t1083 * t842 - t1085 * t841 - t1097 * t840;
t650 = -0.2e1 * t1080 * t843 - 0.2e1 * t1082 * t842 - 0.2e1 * t1084 * t841 - 0.2e1 * t1096 * t840;
t637 = t687 * t1433 + t684 * t865;
t636 = t686 * t1434 + t683 * t865;
t635 = t688 * t1435 + t682 * t865;
t634 = t685 * t1436 - t681 * t865;
t633 = -t1215 * t847 - t1217 * t846 - t1219 * t845 - t1231 * t844;
t632 = t1215 * t843 + t1217 * t842 + t1219 * t841 + t1231 * t840;
t631 = -t1073 * t847 - t1076 * t846 - t1079 * t845 - t1095 * t844;
t630 = t1073 * t843 + t1076 * t842 + t1079 * t841 + t1095 * t840;
t629 = -0.2e1 * t1072 * t847 - 0.2e1 * t1075 * t846 - 0.2e1 * t1078 * t845 - 0.2e1 * t1094 * t844;
t628 = 0.2e1 * t1072 * t843 + 0.2e1 * t1075 * t842 + 0.2e1 * t1078 * t841 + 0.2e1 * t1094 * t840;
t627 = -t1300 * t684 - t637 * t878;
t626 = -t1301 * t684 + t637 * t884;
t625 = -t1302 * t683 - t636 * t876;
t624 = -t1303 * t683 + t636 * t882;
t623 = -t1304 * t682 - t635 * t874;
t622 = -t1305 * t682 + t635 * t880;
t621 = t1306 * t681 - t634 * t866;
t620 = t1307 * t681 + t634 * t868;
t619 = t728 * t1498 + t726 * t1501 + t724 * t1504 + t722 * t1507;
t618 = t727 * t1498 + t725 * t1501 + t723 * t1504 + t721 * t1507;
t617 = t1552 * t721 * t722 + t1554 * t723 * t724 + t1556 * t725 * t726 + t1558 * t727 * t728;
t616 = t658 * t1499 + t657 * t1502 + t656 * t1505 + t653 * t1508;
t615 = (t1179 * t958 + t1185 * t960 + t1191 * t962 + t1197 * t964) * t863;
t614 = (-t1180 * t958 - t1186 * t960 - t1192 * t962 - t1198 * t964) * t863;
t613 = (-t1179 * t959 - t1185 * t961 - t1191 * t963 - t1197 * t965) * t863;
t612 = (t1180 * t959 + t1186 * t961 + t1192 * t963 + t1198 * t965) * t863;
t611 = t658 * t1510 + t657 * t1512 + t656 * t1514 + t653 * t1516;
t610 = t658 * t1511 + t657 * t1513 + t656 * t1515 + t653 * t1517;
t609 = (t1048 * t1057 + t1049 * t1062 + t1050 * t1067 + t1051 * t1088) * t863;
t608 = (-t1048 * t1058 - t1049 * t1063 - t1050 * t1068 - t1051 * t1089) * t863;
t1 = [t1552 * t721 ^ 2 + t1554 * t723 ^ 2 + t1556 * t725 ^ 2 + t1558 * t727 ^ 2, t1312 * t770 + t1323 * t769 + t1334 * t768 + t1365 * t836, (-t1019 * t860 - t1022 * t858 - t1025 * t856 - t1028 * t851) * t1528, (t1018 * t847 + t1021 * t846 + t1024 * t845 + t1027 * t844) * t1528, t1221 * t839 + t1223 * t838 + t1225 * t837 + t1233 * t836, 0.2e1 * t1220 * t839 + 0.2e1 * t1222 * t838 + 0.2e1 * t1224 * t837 + 0.2e1 * t1232 * t836, (-t1040 * t1536 - t1041 * t1537 - t1042 * t1538 + t836 * t939) * t895, (-t1040 * t860 - t1041 * t858 - t1042 * t856 + t1043 * t836) * t895, (t1262 * t836 + t1312 * t1574 + t1323 * t1575 + t1334 * t1576) * t896, t655 * t1517 + t662 * t1515 + t664 * t1513 + t666 * t1511 + (-t1004 * t725 - t1006 * t723 - t1008 * t727 + t721 * t998) * t895 + (-t1019 - t1022 - t1025 - t1028) * t863, t660 * t1517 + t668 * t1515 + t670 * t1513 + t672 * t1511 + (t1003 * t725 + t1005 * t723 + t1007 * t727 - t721 * t996) * t895 + (t721 * t926 + t723 * t924 + t725 * t922 + t727 * t920) * t863, 0, 0, 0, t827; t617, t652, t609, t608, t651, t650, t1585, t1586, t1577, t654 * t1517 + t661 * t1515 + t663 * t1513 + t665 * t1511 + (-t1004 * t726 - t1006 * t724 - t1008 * t728 + t722 * t998) * t895 + (-t1059 * t728 - t1064 * t726 - t1069 * t724 - t1090 * t722) * t863, t659 * t1517 + t667 * t1515 + t669 * t1513 + t671 * t1511 + (t1003 * t726 + t1005 * t724 + t1007 * t728 - t722 * t996) * t895 + (t722 * t926 + t724 * t924 + t726 * t922 + t728 * t920) * t863, 0, 0, 0, 0; t618, t633, t613, t612, t631, t629, t1587, t1588, t1571, t673 * t1517 + t675 * t1515 + t676 * t1513 + t677 * t1511 + (-t1004 * t735 - t1006 * t734 - t1008 * t736 + t730 * t998) * t895 + (-t1059 * t736 - t1064 * t735 - t1069 * t734 - t1090 * t730) * t863, t674 * t1517 + t678 * t1515 + t679 * t1513 + t680 * t1511 + (t1003 * t735 + t1005 * t734 + t1007 * t736 - t730 * t996) * t895 + (t730 * t926 + t734 * t924 + t735 * t922 + t736 * t920) * t863, 0, 0, 0, 0; t610, -t1201 * t687 - t1206 * t686 - t1211 * t688 - t1227 * t685, ((-t1161 + t1388) * t885 + (-t1163 + t1391) * t883 + (-t1165 + t1385) * t881 + (-t1167 + t1395) * t869) * t863, (-t1385 * t875 - t1388 * t879 - t1391 * t877 - t1395 * t867) * t863 + t1173, -t1148 * t845 - t1150 * t847 - t1152 * t846 - t1154 * t844, -0.2e1 * t1281 * t845 - 0.2e1 * t1284 * t847 - 0.2e1 * t1287 * t846 - 0.2e1 * t1290 * t844, t844 * t1159 - t845 * t1158 - t846 * t1157 - t847 * t1156 + (-t1264 * t686 * t876 - t1265 * t688 * t874 - t1266 * t687 * t878 + t1155 * t844) * t895, t844 * t1399 - t845 * t1398 - t846 * t1397 - t847 * t1396 + (-t1518 * t751 * t845 - t1519 * t753 * t847 - t1520 * t752 * t846 + t1394 * t844) * t895, (-t1258 * t681 - t1264 * t683 - t1265 * t682 - t1266 * t684) * t895, t620 * t1517 + t622 * t1515 + t624 * t1513 + t626 * t1511 + (t1124 * t653 - t1131 * t657 - t1133 * t656 - t1135 * t658) * t895 + (-t1200 * t658 - t1205 * t657 - t1210 * t656 - t1226 * t653) * t863, t621 * t1517 + t623 * t1515 + t625 * t1513 + t627 * t1511 + (-t1122 * t653 + t1130 * t657 + t1132 * t656 + t1134 * t658) * t895 + (t653 * t980 + t656 * t974 + t657 * t972 + t658 * t970) * t863, 0, -t848, -t849, 0; t617, t652, t609, t608, t651, t650, t1585, t1586, t1577, t655 * t1516 + t662 * t1514 + t664 * t1512 + t666 * t1510 + (-t721 * t999 + t723 * t993 + t725 * t989 + t727 * t985) * t895 + (t1060 * t727 + t1065 * t725 + t1070 * t723 + t1091 * t721) * t863, t660 * t1516 + t668 * t1514 + t670 * t1512 + t672 * t1510 + (t721 * t997 - t723 * t992 - t725 * t988 - t727 * t984) * t895 + (-t721 * t927 - t723 * t925 - t725 * t923 - t727 * t921) * t863, 0, 0, 0, 0; t1552 * t722 ^ 2 + t1554 * t724 ^ 2 + t1556 * t726 ^ 2 + t1558 * t728 ^ 2, t1313 * t770 + t1324 * t769 + t1335 * t768 + t1365 * t832, (t1017 * t860 + t1020 * t858 + t1023 * t856 + t1026 * t851) * t1528, (-t1058 * t728 * t843 - t1063 * t726 * t842 - t1068 * t724 * t841 - t1089 * t722 * t840) * t1528, t1221 * t835 + t1223 * t834 + t1225 * t833 + t1233 * t832, 0.2e1 * t1220 * t835 + 0.2e1 * t1222 * t834 + 0.2e1 * t1224 * t833 + 0.2e1 * t1232 * t832, (t1168 * t1536 + t1169 * t1537 + t1170 * t1538 + t832 * t939) * t895, (t1043 * t832 + t1168 * t860 + t1169 * t858 + t1170 * t856) * t895, (t1262 * t832 + t1313 * t1574 + t1324 * t1575 + t1335 * t1576) * t896, t654 * t1516 + t661 * t1514 + t663 * t1512 + t665 * t1510 + (-t722 * t999 + t724 * t993 + t726 * t989 + t728 * t985) * t895 + (t1017 + t1020 + t1023 + t1026) * t863, t659 * t1516 + t667 * t1514 + t669 * t1512 + t671 * t1510 + (t722 * t997 - t724 * t992 - t726 * t988 - t728 * t984) * t895 + (-t722 * t927 - t724 * t925 - t726 * t923 - t728 * t921) * t863, 0, 0, 0, t827; t619, t632, t615, t614, t630, t628, t1564, t1563, t1562, t673 * t1516 + t675 * t1514 + t676 * t1512 + t677 * t1510 + (-t730 * t999 + t734 * t993 + t735 * t989 + t736 * t985) * t895 + (t1060 * t736 + t1065 * t735 + t1070 * t734 + t1091 * t730) * t863, t674 * t1516 + t678 * t1514 + t679 * t1512 + t680 * t1510 + (t730 * t997 - t734 * t992 - t735 * t988 - t736 * t984) * t895 + (-t730 * t927 - t734 * t925 - t735 * t923 - t736 * t921) * t863, 0, 0, 0, 0; t611, t1203 * t687 + t1208 * t686 + t1213 * t688 + t1229 * t685, (t1383 * t724 + t1386 * t728 + t1389 * t726 + t1392 * t722) * t863 + t1174, ((-t1203 * t658 - t687 * t1510) * t879 + (-t1208 * t657 - t686 * t1512) * t877 + (-t1213 * t656 - t688 * t1514) * t875 + (-t1229 * t653 - t685 * t1516) * t867) * t863, t1148 * t841 + t1150 * t843 + t1152 * t842 + t1154 * t840, 0.2e1 * t1281 * t841 + 0.2e1 * t1284 * t843 + 0.2e1 * t1287 * t842 + 0.2e1 * t1290 * t840, -t840 * t1159 + t841 * t1158 + t842 * t1157 + t843 * t1156 + (t1149 * t841 + t1151 * t843 + t1153 * t842 - t1155 * t840) * t895, -t840 * t1399 + t841 * t1398 + t842 * t1397 + t843 * t1396 + (t1485 * t687 * t843 + t1487 * t686 * t842 + t1489 * t688 * t841 - t1394 * t840) * t895, (t1243 * t684 + t1246 * t683 + t1249 * t682 + t1259 * t681) * t895, t620 * t1516 + t622 * t1514 + t624 * t1512 + t626 * t1510 + (t1106 * t658 + t1111 * t657 + t1116 * t656 - t1125 * t653) * t895 + (t1202 * t658 + t1207 * t657 + t1212 * t656 + t1228 * t653) * t863, t621 * t1516 + t623 * t1514 + t625 * t1512 + t627 * t1510 + (-t1105 * t658 - t1110 * t657 - t1115 * t656 + t1123 * t653) * t895 + (-t653 * t981 - t656 * t975 - t657 * t973 - t658 * t971) * t863, 0, t849, -t848, 0; t618, t633, t613, t612, t631, t629, t1587, t1588, t1571, t655 * t1508 + t662 * t1505 + t664 * t1502 + t666 * t1499 + (t1099 * t727 + t1101 * t725 + t1103 * t723 + t1121 * t721) * t895 + (t1214 * t727 + t1216 * t725 + t1218 * t723 + t1230 * t721) * t863, t660 * t1508 + t668 * t1505 + t670 * t1502 + t672 * t1499 + (-t1098 * t727 - t1100 * t725 - t1102 * t723 - t1120 * t721) * t895 + (-t721 * t983 - t723 * t978 - t725 * t977 - t727 * t976) * t863, 0, 0, 0, 0; t619, t632, t615, t614, t630, t628, t1564, t1563, t1562, t654 * t1508 + t661 * t1505 + t663 * t1502 + t665 * t1499 + (t1099 * t728 + t1101 * t726 + t1103 * t724 + t1121 * t722) * t895 + (t1214 * t728 + t1216 * t726 + t1218 * t724 + t1230 * t722) * t863, t659 * t1508 + t667 * t1505 + t669 * t1502 + t671 * t1499 + (-t1098 * t728 - t1100 * t726 - t1102 * t724 - t1120 * t722) * t895 + (-t722 * t983 - t724 * t978 - t726 * t977 - t728 * t976) * t863, 0, 0, 0, 0; t1552 * t730 ^ 2 + t1554 * t734 ^ 2 + t1556 * t735 ^ 2 + t1558 * t736 ^ 2, t1446 * t767 + t1452 * t766 + t1458 * t765 + t1464 * t761, (t1144 * t860 + t1145 * t858 + t1146 * t856 + t1147 * t851) * t1528, (-t1180 * t736 * t773 - t1186 * t735 * t772 - t1192 * t734 * t771 - t1198 * t730 * t763) * t1528, t1310 * t767 + t1321 * t766 + t1332 * t765 + t1343 * t761, 0.2e1 * t1309 * t767 + 0.2e1 * t1320 * t766 + 0.2e1 * t1331 * t765 + 0.2e1 * t1342 * t761, (t1308 * t1544 + t1319 * t1545 + t1330 * t1546 + t1341 * t1547) * t1527, (t1447 * t1544 + t1453 * t1545 + t1459 * t1546 + t1465 * t1547) * t1527, (t1446 * t756 ^ 2 + t1452 * t755 ^ 2 + t1458 * t754 ^ 2 + t1464 * t747 ^ 2) * t896, t673 * t1508 + t675 * t1505 + t676 * t1502 + t677 * t1499 + (t1099 * t736 + t1101 * t735 + t1103 * t734 + t1121 * t730) * t895 + (t1144 + t1145 + t1146 + t1147) * t863, t674 * t1508 + t678 * t1505 + t679 * t1502 + t680 * t1499 + (-t1098 * t736 - t1100 * t735 - t1102 * t734 - t1120 * t730) * t895 + (-t730 * t983 - t734 * t978 - t735 * t977 - t736 * t976) * t863, 0, 0, 0, 1; t616, t1357 * t687 + t1359 * t686 + t1361 * t688 + t1364 * t685, (t1383 * t734 + t1386 * t736 + t1389 * t735 + t1392 * t730) * t863 + t1172, ((-t1357 * t658 - t687 * t1499) * t879 + (-t1359 * t657 - t686 * t1502) * t877 + (-t1361 * t656 - t688 * t1505) * t875 + (-t1364 * t653 - t685 * t1508) * t867) * t863, t1280 * t771 + t1283 * t773 + t1286 * t772 + t1289 * t763, 0.2e1 * t1384 * t771 + 0.2e1 * t1387 * t773 + 0.2e1 * t1390 * t772 + 0.2e1 * t1393 * t763, -t763 * t1294 + t771 * t1293 + t772 * t1292 + t773 * t1291 + (t1279 * t754 + t1282 * t756 + t1285 * t755 + t1288 * t747) * t895, -t763 * t1551 * t681 + t771 * t1553 * t682 + t772 * t1555 * t683 + t773 * t1557 * t684 + (t1518 * t754 + t1519 * t756 + t1520 * t755 + t1521 * t747) * t895, (t1367 * t684 + t1368 * t683 + t1369 * t682 - t1376 * t681) * t895, t620 * t1508 + t622 * t1505 + t624 * t1502 + t626 * t1499 + (t1236 * t658 + t1239 * t657 + t1242 * t656 + t1257 * t653) * t895 + (t1356 * t658 + t1358 * t657 + t1360 * t656 + t1363 * t653) * t863, t621 * t1508 + t623 * t1505 + t625 * t1502 + t627 * t1499 + (-t1235 * t658 - t1238 * t657 - t1241 * t656 - t1256 * t653) * t895 + (-t1071 * t658 - t1074 * t657 - t1077 * t656 - t1093 * t653) * t863, 0, 0, 0, 0; t610, -t1269 * t846 - t1272 * t845 - t1275 * t844 - t1278 * t847, ((t727 * t1497 - t658 * t847) * t776 * t1314 + (t725 * t1494 - t657 * t846) * t775 * t1325 + (t723 * t1495 - t656 * t845) * t774 * t1336 + (t721 * t1496 - t653 * t844) * t764 * t1345) * t863, (-t1018 * t1044 - t1021 * t1045 - t1024 * t1046 - t1027 * t1047) * t863 + t1173, -t1044 * t1081 - t1045 * t1083 - t1046 * t1085 - t1047 * t1097, -0.2e1 * t1044 * t1080 - 0.2e1 * t1045 * t1082 - 0.2e1 * t1046 * t1084 - 0.2e1 * t1047 * t1096, (-t1536 * t940 - t1537 * t942 - t1538 * t944 + t1539 * t946) * t895, (t851 * t946 - t856 * t944 - t858 * t942 - t860 * t940) * t895, (-t1252 * t847 - t1253 * t846 - t1254 * t845 - t1263 * t844) * t896, t653 * t655 + t656 * t662 + t657 * t664 + t658 * t666 + (-t1002 * t721 + t723 * t995 + t725 * t991 + t727 * t987) * t895 + (t1136 * t725 + t1138 * t723 + t1140 * t721 + t1142 * t727) * t863, t653 * t660 + t656 * t668 + t657 * t670 + t658 * t672 + (t1001 * t721 - t723 * t994 - t725 * t990 - t727 * t986) * t895 + (-t721 * t934 - t723 * t933 - t725 * t932 - t727 * t935) * t863, 0, -t848, -t849, 0; t611, t1137 + t1139 + t1141 + t1143, (t1010 * t726 + t1012 * t724 + t1014 * t722 + t1016 * t728) * t863 + t1174, ((-t728 * t1497 - t658 * t843) * t776 * t1315 + (-t726 * t1494 - t657 * t842) * t775 * t1326 + (-t724 * t1495 - t656 * t841) * t774 * t1337 + (-t722 * t1496 - t653 * t840) * t764 * t1346) * t863, t1137 * t854 + t1139 * t853 + t1141 * t850 + t1143 * t855, 0.2e1 * t1044 * t1220 * t843 + 0.2e1 * t1045 * t1222 * t842 + 0.2e1 * t1046 * t1224 * t841 + 0.2e1 * t1047 * t1232 * t840, (t1536 * t941 + t1537 * t943 + t1538 * t945 - t1539 * t947) * t895, (-t851 * t947 + t856 * t945 + t858 * t943 + t860 * t941) * t895, (t1252 * t843 + t1253 * t842 + t1254 * t841 + t1263 * t840) * t896, t653 * t654 + t656 * t661 + t657 * t663 + t658 * t665 + (-t1002 * t722 + t724 * t995 + t726 * t991 + t728 * t987) * t895 + (t1136 * t726 + t1138 * t724 + t1140 * t722 + t1142 * t728) * t863, t653 * t659 + t656 * t667 + t657 * t669 + t658 * t671 + (t1001 * t722 - t724 * t994 - t726 * t990 - t728 * t986) * t895 + (-t722 * t934 - t724 * t933 - t726 * t932 - t728 * t935) * t863, 0, t849, -t848, 0; t616, t1044 * t1215 + t1045 * t1217 + t1046 * t1219 + t1047 * t1231, (t1010 * t735 + t1012 * t734 + t1014 * t730 + t1016 * t736) * t863 + t1172, ((-t1497 * t1500 - t658 * t773) * t1315 + (-t1494 * t1503 - t657 * t772) * t1326 + (-t1495 * t1506 - t656 * t771) * t1337 + (-t1496 * t1509 - t653 * t763) * t1346) * t863, t1044 * t1073 + t1045 * t1076 + t1046 * t1079 + t1047 * t1095, 0.2e1 * t1044 * t1072 + 0.2e1 * t1045 * t1075 + 0.2e1 * t1046 * t1078 + 0.2e1 * t1047 * t1094, (t1536 * t901 + t1537 * t902 + t1538 * t903 + t1539 * t904) * t895, (t851 * t904 + t856 * t903 + t858 * t902 + t860 * t901) * t895, (t1044 * t1245 + t1045 * t1248 + t1046 * t1251 - t1366 * t1573) * t896, t653 * t673 + t656 * t675 + t657 * t676 + t658 * t677 + (-t1002 * t730 + t734 * t995 + t735 * t991 + t736 * t987) * t895 + (t1136 * t735 + t1138 * t734 + t1140 * t730 + t1142 * t736) * t863, t653 * t674 + t656 * t678 + t657 * t679 + t658 * t680 + (t1001 * t730 - t734 * t994 - t735 * t990 - t736 * t986) * t895 + (-t730 * t934 - t734 * t933 - t735 * t932 - t736 * t935) * t863, 0, 0, 0, 0; t653 ^ 2 + t656 ^ 2 + t657 ^ 2 + t658 ^ 2, t1268 * t686 + t1271 * t688 + t1274 * t685 + t1277 * t687, (t869 * t1532 + t883 * t1533 + t885 * t1534 + t881 * t1535) * t863, (-t867 * t1532 - t877 * t1533 - t879 * t1534 - t875 * t1535) * t863, t1044 * t1150 + t1045 * t1152 + t1046 * t1148 + t1047 * t1154, 0.2e1 * t1044 * t1284 + 0.2e1 * t1045 * t1287 + 0.2e1 * t1046 * t1281 + 0.2e1 * t1047 * t1290, -t1047 * t1159 + t1046 * t1158 + t1045 * t1157 + t1044 * t1156 + (-t1261 * t685 * t866 + t1044 * t1151 + t1045 * t1153 + t1046 * t1149) * t895, -t681 * t1381 + t682 * t1380 + t683 * t1379 + t684 * t1382 + (t1371 * t687 + t1373 * t686 + t1375 * t688 - t1378 * t685) * t895, (t1244 * t684 + t1247 * t683 + t1250 * t682 + t1261 * t681) * t895, t653 * t620 + t656 * t622 + t657 * t624 + t658 * t626 + (t1108 * t658 + t1113 * t657 + t1118 * t656 - t1129 * t653) * t895 + (t1267 * t657 + t1270 * t656 + t1273 * t653 + t1276 * t658) * t863, t653 * t621 + t656 * t623 + t657 * t625 + t658 * t627 + (-t1107 * t658 - t1112 * t657 - t1117 * t656 + t1128 * t653) * t895 + (-t1009 * t657 - t1011 * t656 - t1013 * t653 - t1015 * t658) * t863, 1, 0, 0, 0;];
tau_reg  = t1;
