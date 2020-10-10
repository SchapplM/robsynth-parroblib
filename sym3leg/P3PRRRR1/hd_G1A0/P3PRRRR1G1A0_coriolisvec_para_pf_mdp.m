% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRRR1G1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR1G1A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR1G1A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR1G1A0_coriolisvec_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:18
% EndTime: 2020-03-09 20:34:20
% DurationCPUTime: 2.55s
% Computational Cost: add. (1923->234), mult. (5706->526), div. (1878->17), fcn. (5343->18), ass. (0->225)
t1501 = legFrame(3,3);
t1468 = sin(t1501);
t1471 = cos(t1501);
t1516 = xDP(2);
t1517 = xDP(1);
t1462 = t1468 * t1517 - t1471 * t1516;
t1510 = cos(qJ(3,3));
t1504 = sin(qJ(3,3));
t1511 = cos(qJ(2,3));
t1618 = t1504 * t1511;
t1441 = (-t1468 * t1516 - t1471 * t1517) * t1510 - t1462 * t1618;
t1505 = sin(qJ(2,3));
t1475 = 0.1e1 / t1505;
t1476 = 0.1e1 / t1505 ^ 2;
t1486 = 0.1e1 / t1510;
t1529 = t1510 ^ 2;
t1487 = 0.1e1 / t1529;
t1488 = t1486 * t1487;
t1636 = t1475 * t1511;
t1591 = t1441 * t1636;
t1619 = t1504 * t1505;
t1518 = 0.1e1 / pkin(2);
t1666 = t1518 ^ 2;
t1405 = (-(-t1462 * t1619 + t1591) * t1476 * t1441 - (-t1441 * t1504 + t1462 * t1511) * t1475 * t1462) * t1488 * t1486 * t1666;
t1489 = 0.1e1 / t1529 ^ 2;
t1438 = t1441 ^ 2;
t1519 = 0.1e1 / pkin(2) ^ 2;
t1596 = t1438 * t1476 * t1519;
t1432 = t1489 * t1596;
t1459 = t1462 ^ 2;
t1642 = t1459 * t1487;
t1588 = t1519 * t1642;
t1426 = t1432 + t1588;
t1665 = -0.2e1 * t1519;
t1546 = t1462 * t1591 * t1665;
t1613 = t1510 * t1511;
t1669 = -t1426 * t1505 * t1510 + t1504 * t1488 * t1546 + t1405 * t1613;
t1502 = legFrame(2,3);
t1469 = sin(t1502);
t1472 = cos(t1502);
t1463 = t1469 * t1517 - t1472 * t1516;
t1512 = cos(qJ(3,2));
t1506 = sin(qJ(3,2));
t1513 = cos(qJ(2,2));
t1616 = t1506 * t1513;
t1442 = (-t1469 * t1516 - t1472 * t1517) * t1512 - t1463 * t1616;
t1507 = sin(qJ(2,2));
t1479 = 0.1e1 / t1507;
t1480 = 0.1e1 / t1507 ^ 2;
t1491 = 0.1e1 / t1512;
t1533 = t1512 ^ 2;
t1492 = 0.1e1 / t1533;
t1493 = t1491 * t1492;
t1633 = t1479 * t1513;
t1590 = t1442 * t1633;
t1617 = t1506 * t1507;
t1406 = (-(-t1463 * t1617 + t1590) * t1480 * t1442 - (-t1442 * t1506 + t1463 * t1513) * t1479 * t1463) * t1493 * t1491 * t1666;
t1494 = 0.1e1 / t1533 ^ 2;
t1439 = t1442 ^ 2;
t1594 = t1439 * t1480 * t1519;
t1433 = t1494 * t1594;
t1460 = t1463 ^ 2;
t1640 = t1460 * t1492;
t1586 = t1519 * t1640;
t1427 = t1433 + t1586;
t1545 = t1463 * t1590 * t1665;
t1612 = t1512 * t1513;
t1668 = -t1427 * t1507 * t1512 + t1506 * t1493 * t1545 + t1406 * t1612;
t1503 = legFrame(1,3);
t1470 = sin(t1503);
t1473 = cos(t1503);
t1464 = t1470 * t1517 - t1473 * t1516;
t1514 = cos(qJ(3,1));
t1508 = sin(qJ(3,1));
t1515 = cos(qJ(2,1));
t1614 = t1508 * t1515;
t1443 = (-t1470 * t1516 - t1473 * t1517) * t1514 - t1464 * t1614;
t1509 = sin(qJ(2,1));
t1483 = 0.1e1 / t1509;
t1484 = 0.1e1 / t1509 ^ 2;
t1496 = 0.1e1 / t1514;
t1537 = t1514 ^ 2;
t1497 = 0.1e1 / t1537;
t1498 = t1496 * t1497;
t1630 = t1483 * t1515;
t1589 = t1443 * t1630;
t1615 = t1508 * t1509;
t1407 = (-(-t1464 * t1615 + t1589) * t1484 * t1443 - (-t1443 * t1508 + t1464 * t1515) * t1483 * t1464) * t1498 * t1496 * t1666;
t1499 = 0.1e1 / t1537 ^ 2;
t1440 = t1443 ^ 2;
t1592 = t1440 * t1484 * t1519;
t1434 = t1499 * t1592;
t1461 = t1464 ^ 2;
t1638 = t1461 * t1497;
t1584 = t1519 * t1638;
t1428 = t1434 + t1584;
t1544 = t1464 * t1589 * t1665;
t1611 = t1514 * t1515;
t1667 = -t1428 * t1509 * t1514 + t1508 * t1498 * t1544 + t1407 * t1611;
t1664 = MDP(2) * t1518;
t1663 = MDP(8) * t1518;
t1520 = t1518 * t1519;
t1662 = MDP(9) * t1520;
t1661 = t1405 * t1475;
t1660 = t1405 * t1511;
t1659 = t1406 * t1479;
t1658 = t1406 * t1513;
t1657 = t1407 * t1483;
t1656 = t1407 * t1515;
t1477 = t1475 * t1476;
t1643 = t1459 * t1475;
t1420 = (t1438 * t1477 + t1643) * t1488 * t1518;
t1655 = t1420 * t1486;
t1654 = t1420 * t1487;
t1481 = t1479 * t1480;
t1641 = t1460 * t1479;
t1421 = (t1439 * t1481 + t1641) * t1493 * t1518;
t1653 = t1421 * t1491;
t1652 = t1421 * t1492;
t1485 = t1483 * t1484;
t1639 = t1461 * t1483;
t1422 = (t1440 * t1485 + t1639) * t1498 * t1518;
t1651 = t1422 * t1496;
t1650 = t1422 * t1497;
t1649 = t1438 * t1511;
t1648 = t1439 * t1513;
t1647 = t1440 * t1515;
t1646 = t1441 * t1462;
t1645 = t1442 * t1463;
t1644 = t1443 * t1464;
t1637 = t1475 * t1486;
t1490 = t1486 * t1489;
t1635 = t1476 * t1490;
t1634 = t1479 * t1491;
t1495 = t1491 * t1494;
t1632 = t1480 * t1495;
t1631 = t1483 * t1496;
t1500 = t1496 * t1499;
t1629 = t1484 * t1500;
t1628 = t1486 * t1405;
t1626 = t1489 * t1504;
t1625 = t1491 * t1406;
t1623 = t1494 * t1506;
t1622 = t1496 * t1407;
t1620 = t1499 * t1508;
t1411 = t1426 * t1619 + t1487 * t1546;
t1412 = t1427 * t1617 + t1492 * t1545;
t1413 = t1428 * t1615 + t1497 * t1544;
t1610 = 0.2e1 * t1520;
t1609 = 0.2e1 * t1646;
t1608 = 0.2e1 * t1645;
t1607 = 0.2e1 * t1644;
t1606 = t1487 * t1661;
t1605 = t1504 * t1628;
t1604 = t1492 * t1659;
t1603 = t1506 * t1625;
t1602 = t1497 * t1657;
t1601 = t1508 * t1622;
t1600 = t1420 * t1637;
t1599 = t1421 * t1634;
t1598 = t1422 * t1631;
t1597 = t1438 * t1635;
t1595 = t1439 * t1632;
t1593 = t1440 * t1629;
t1587 = t1459 * t1626;
t1585 = t1460 * t1623;
t1583 = t1461 * t1620;
t1582 = t1487 * t1636;
t1581 = t1476 * t1626;
t1580 = t1492 * t1633;
t1579 = t1480 * t1623;
t1578 = t1497 * t1630;
t1577 = t1484 * t1620;
t1576 = 0.2e1 * t1504 * t1661;
t1575 = 0.2e1 * t1506 * t1659;
t1574 = 0.2e1 * t1508 * t1657;
t1474 = t1504 ^ 2;
t1573 = t1474 * t1606;
t1478 = t1506 ^ 2;
t1572 = t1478 * t1604;
t1482 = t1508 ^ 2;
t1571 = t1482 * t1602;
t1570 = t1420 * t1582;
t1569 = t1421 * t1580;
t1568 = t1422 * t1578;
t1567 = t1477 * t1490 * t1649;
t1566 = t1481 * t1495 * t1648;
t1565 = t1485 * t1500 * t1647;
t1564 = t1459 * t1474 * t1488 * t1505;
t1563 = t1460 * t1478 * t1493 * t1507;
t1562 = t1461 * t1482 * t1498 * t1509;
t1561 = t1628 * t1636;
t1560 = t1504 * t1582;
t1559 = t1625 * t1633;
t1558 = t1506 * t1580;
t1557 = t1622 * t1630;
t1556 = t1508 * t1578;
t1555 = (-t1519 * t1564 + t1669) * t1637;
t1554 = (-t1519 * t1563 + t1668) * t1634;
t1553 = (-t1519 * t1562 + t1667) * t1631;
t1552 = ((-t1505 * t1588 - t1660) * t1504 + t1411) * t1637;
t1551 = ((-t1507 * t1586 - t1658) * t1506 + t1412) * t1634;
t1550 = ((-t1509 * t1584 - t1656) * t1508 + t1413) * t1631;
t1549 = (-0.1e1 + 0.2e1 * t1529) * t1635 * t1646;
t1548 = (-0.1e1 + 0.2e1 * t1533) * t1632 * t1645;
t1547 = (-0.1e1 + 0.2e1 * t1537) * t1629 * t1644;
t1543 = (t1474 * t1490 + t1488) * t1643;
t1542 = (t1478 * t1495 + t1493) * t1641;
t1541 = (t1482 * t1500 + t1498) * t1639;
t1458 = t1470 * t1611 - t1473 * t1508;
t1457 = t1469 * t1612 - t1472 * t1506;
t1456 = t1468 * t1613 - t1471 * t1504;
t1455 = -t1470 * t1514 + t1473 * t1614;
t1454 = t1470 * t1508 + t1473 * t1611;
t1453 = -t1469 * t1512 + t1472 * t1616;
t1452 = t1469 * t1506 + t1472 * t1612;
t1451 = -t1468 * t1510 + t1471 * t1618;
t1450 = t1468 * t1504 + t1471 * t1613;
t1449 = -t1470 * t1614 - t1473 * t1514;
t1448 = -t1469 * t1616 - t1472 * t1512;
t1447 = -t1468 * t1618 - t1471 * t1510;
t1416 = -0.2e1 * t1497 * t1592 + t1434;
t1415 = -0.2e1 * t1492 * t1594 + t1433;
t1414 = -0.2e1 * t1487 * t1596 + t1432;
t1 = [(t1450 * t1600 + t1452 * t1599 + t1454 * t1598) * MDP(1) + (t1447 * t1606 + t1448 * t1604 + t1449 * t1602) * t1664 + (t1450 * t1561 + t1452 * t1559 + t1454 * t1557 + (-t1450 * t1597 - t1452 * t1595 - t1454 * t1593) * t1519 + (t1447 * t1570 + t1448 * t1569 + t1449 * t1568) * t1518) * MDP(3) + (-t1450 * t1628 - t1452 * t1625 - t1454 * t1622 + (-t1450 * t1567 - t1452 * t1566 - t1454 * t1565) * t1519 + (-t1447 * t1654 - t1448 * t1652 - t1449 * t1650) * t1518) * MDP(4) + ((t1447 * t1573 + t1448 * t1572 + t1449 * t1571) * t1518 + ((-t1440 * t1470 + t1449 * t1607) * t1577 + (-t1439 * t1469 + t1448 * t1608) * t1579 + (-t1438 * t1468 + t1447 * t1609) * t1581) * t1520) * MDP(5) + ((t1447 * t1549 + t1448 * t1548 + t1449 * t1547) * t1610 + ((t1416 * t1470 + t1449 * t1574) * t1496 + (t1415 * t1469 + t1448 * t1575) * t1491 + (t1414 * t1468 + t1447 * t1576) * t1486) * t1518) * MDP(6) + ((t1468 * t1605 + t1469 * t1603 + t1470 * t1601) * t1518 + (t1447 * t1543 + t1448 * t1542 + t1449 * t1541) * t1520) * MDP(7) + (t1405 * t1468 + t1406 * t1469 + t1407 * t1470) * t1663 + (t1468 * t1587 + t1469 * t1585 + t1470 * t1583) * t1662 + (t1454 * t1553 + t1452 * t1554 + t1450 * t1555 + ((t1449 * t1630 - t1470 * t1615) * t1651 + (t1448 * t1633 - t1469 * t1617) * t1653 + (t1447 * t1636 - t1468 * t1619) * t1655) * t1518) * MDP(10) + (t1454 * t1550 + t1452 * t1551 + t1450 * t1552 + ((-t1449 * t1556 - t1470 * t1509) * t1422 + (-t1448 * t1558 - t1469 * t1507) * t1421 + (-t1447 * t1560 - t1468 * t1505) * t1420) * t1518) * MDP(11); (t1456 * t1600 + t1457 * t1599 + t1458 * t1598) * MDP(1) + (t1451 * t1606 + t1453 * t1604 + t1455 * t1602) * t1664 + (t1456 * t1561 + t1457 * t1559 + t1458 * t1557 + (-t1456 * t1597 - t1457 * t1595 - t1458 * t1593) * t1519 + (t1451 * t1570 + t1453 * t1569 + t1455 * t1568) * t1518) * MDP(3) + (-t1456 * t1628 - t1457 * t1625 - t1458 * t1622 + (-t1456 * t1567 - t1457 * t1566 - t1458 * t1565) * t1519 + (-t1451 * t1654 - t1453 * t1652 - t1455 * t1650) * t1518) * MDP(4) + ((t1451 * t1573 + t1453 * t1572 + t1455 * t1571) * t1518 + ((t1440 * t1473 + t1455 * t1607) * t1577 + (t1439 * t1472 + t1453 * t1608) * t1579 + (t1438 * t1471 + t1451 * t1609) * t1581) * t1520) * MDP(5) + ((t1451 * t1549 + t1453 * t1548 + t1455 * t1547) * t1610 + ((-t1416 * t1473 + t1455 * t1574) * t1496 + (-t1415 * t1472 + t1453 * t1575) * t1491 + (-t1414 * t1471 + t1451 * t1576) * t1486) * t1518) * MDP(6) + ((-t1471 * t1605 - t1472 * t1603 - t1473 * t1601) * t1518 + (t1451 * t1543 + t1453 * t1542 + t1455 * t1541) * t1520) * MDP(7) + (-t1405 * t1471 - t1406 * t1472 - t1407 * t1473) * t1663 + (-t1471 * t1587 - t1472 * t1585 - t1473 * t1583) * t1662 + (t1458 * t1553 + t1457 * t1554 + t1456 * t1555 + ((t1455 * t1630 + t1473 * t1615) * t1651 + (t1453 * t1633 + t1472 * t1617) * t1653 + (t1451 * t1636 + t1471 * t1619) * t1655) * t1518) * MDP(10) + (t1458 * t1550 + t1457 * t1551 + t1456 * t1552 + ((-t1455 * t1556 + t1473 * t1509) * t1422 + (-t1453 * t1558 + t1472 * t1507) * t1421 + (-t1451 * t1560 + t1471 * t1505) * t1420) * t1518) * MDP(11); (t1422 + t1421 + t1420) * MDP(1) + (t1656 + t1658 + t1660) * MDP(3) + (-t1405 * t1505 - t1406 * t1507 - t1407 * t1509) * MDP(4) + (t1667 + t1668 + t1669) * MDP(10) + (-t1405 * t1618 - t1406 * t1616 - t1407 * t1614 + t1411 + t1412 + t1413) * MDP(11) + ((-t1438 * t1475 * t1489 - t1439 * t1479 * t1494 - t1440 * t1483 * t1499) * MDP(3) + (-t1476 * t1489 * t1649 - t1480 * t1494 * t1648 - t1484 * t1499 * t1647) * MDP(4) + (-t1562 - t1563 - t1564) * MDP(10) + (-t1615 * t1638 - t1617 * t1640 - t1619 * t1642) * MDP(11)) * t1519;];
taucX  = t1;
