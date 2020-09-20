% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPRR1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRR1G2A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RPRR1G2A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2A0_inertia_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3RPRR1G2A0_inertia_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:25:12
% EndTime: 2020-03-09 21:25:14
% DurationCPUTime: 2.39s
% Computational Cost: add. (5739->294), mult. (3564->461), div. (669->7), fcn. (2937->62), ass. (0->219)
t1570 = 0.1e1 / pkin(3);
t1631 = t1570 / 0.2e1;
t1571 = pkin(1) ^ 2;
t1671 = MDP(4) * t1571 + MDP(1);
t1680 = t1671 / 0.4e1;
t1679 = MDP(5) / 0.2e1;
t1678 = MDP(6) / 0.2e1;
t1677 = MDP(7) / 0.2e1;
t1630 = t1570 / 0.4e1;
t1620 = 0.2e1 * MDP(7);
t1676 = t1620 / 0.2e1;
t1564 = sin(qJ(3,1));
t1623 = pkin(7) + qJ(3,1);
t1513 = pkin(1) * sin(t1623) + t1564 * pkin(2);
t1510 = 0.1e1 / t1513 ^ 2;
t1543 = cos(qJ(1,1) + t1623);
t1635 = t1510 * t1543;
t1675 = t1635 / 0.2e1;
t1563 = sin(qJ(3,2));
t1622 = pkin(7) + qJ(3,2);
t1512 = pkin(1) * sin(t1622) + t1563 * pkin(2);
t1508 = 0.1e1 / t1512 ^ 2;
t1542 = cos(qJ(1,2) + t1622);
t1637 = t1508 * t1542;
t1674 = t1637 / 0.2e1;
t1562 = sin(qJ(3,3));
t1621 = pkin(7) + qJ(3,3);
t1511 = pkin(1) * sin(t1621) + t1562 * pkin(2);
t1506 = 0.1e1 / t1511 ^ 2;
t1541 = cos(qJ(1,3) + t1621);
t1639 = t1506 * t1541;
t1673 = t1639 / 0.2e1;
t1672 = MDP(5) * t1631;
t1557 = sin(pkin(7));
t1670 = pkin(1) * t1557;
t1559 = legFrame(3,2);
t1627 = qJ(1,3) + pkin(7);
t1535 = t1559 + t1627;
t1526 = qJ(3,3) + t1535;
t1536 = -t1559 + t1627;
t1527 = qJ(3,3) + t1536;
t1490 = sin(t1526) + sin(t1527);
t1545 = qJ(1,3) + t1559;
t1546 = qJ(1,3) - t1559;
t1454 = t1490 * pkin(3) + (sin(t1535) + sin(t1536)) * pkin(2) + (sin(t1545) + sin(t1546)) * pkin(1);
t1505 = 0.1e1 / t1511;
t1667 = t1454 * t1505;
t1560 = legFrame(2,2);
t1628 = qJ(1,2) + pkin(7);
t1537 = t1560 + t1628;
t1528 = qJ(3,2) + t1537;
t1538 = -t1560 + t1628;
t1529 = qJ(3,2) + t1538;
t1491 = sin(t1528) + sin(t1529);
t1547 = qJ(1,2) + t1560;
t1548 = qJ(1,2) - t1560;
t1455 = t1491 * pkin(3) + (sin(t1537) + sin(t1538)) * pkin(2) + (sin(t1547) + sin(t1548)) * pkin(1);
t1507 = 0.1e1 / t1512;
t1666 = t1455 * t1507;
t1561 = legFrame(1,2);
t1629 = qJ(1,1) + pkin(7);
t1539 = t1561 + t1629;
t1530 = qJ(3,1) + t1539;
t1540 = -t1561 + t1629;
t1531 = qJ(3,1) + t1540;
t1492 = sin(t1530) + sin(t1531);
t1549 = qJ(1,1) + t1561;
t1550 = qJ(1,1) - t1561;
t1456 = t1492 * pkin(3) + (sin(t1539) + sin(t1540)) * pkin(2) + (sin(t1549) + sin(t1550)) * pkin(1);
t1509 = 0.1e1 / t1513;
t1665 = t1456 * t1509;
t1493 = -cos(t1527) + cos(t1526);
t1457 = -t1493 * pkin(3) + (cos(t1536) - cos(t1535)) * pkin(2) + (cos(t1546) - cos(t1545)) * pkin(1);
t1664 = t1457 * t1505;
t1494 = -cos(t1529) + cos(t1528);
t1458 = -t1494 * pkin(3) + (cos(t1538) - cos(t1537)) * pkin(2) + (cos(t1548) - cos(t1547)) * pkin(1);
t1663 = t1458 * t1507;
t1495 = -cos(t1531) + cos(t1530);
t1459 = -t1495 * pkin(3) + (cos(t1540) - cos(t1539)) * pkin(2) + (cos(t1550) - cos(t1549)) * pkin(1);
t1662 = t1459 * t1509;
t1487 = -pkin(3) * t1541 - pkin(2) * cos(t1627) - cos(qJ(1,3)) * pkin(1);
t1658 = t1487 * t1505;
t1475 = t1570 * t1658;
t1496 = t1541 * t1505;
t1464 = t1496 + t1475;
t1661 = t1464 * t1505;
t1488 = -pkin(3) * t1542 - pkin(2) * cos(t1628) - cos(qJ(1,2)) * pkin(1);
t1657 = t1488 * t1507;
t1476 = t1570 * t1657;
t1497 = t1542 * t1507;
t1466 = t1497 + t1476;
t1660 = t1466 * t1507;
t1489 = -pkin(3) * t1543 - pkin(2) * cos(t1629) - cos(qJ(1,1)) * pkin(1);
t1656 = t1489 * t1509;
t1477 = t1570 * t1656;
t1498 = t1543 * t1509;
t1468 = t1498 + t1477;
t1659 = t1468 * t1509;
t1655 = t1490 * t1505;
t1654 = t1491 * t1507;
t1653 = t1492 * t1509;
t1652 = t1493 * t1505;
t1651 = t1494 * t1507;
t1650 = t1495 * t1509;
t1558 = cos(pkin(7));
t1544 = pkin(1) * t1558 + pkin(2);
t1565 = cos(qJ(3,3));
t1634 = t1557 * t1562;
t1499 = -pkin(1) * t1634 + t1544 * t1565;
t1649 = t1499 * t1506;
t1566 = cos(qJ(3,2));
t1633 = t1557 * t1563;
t1500 = -pkin(1) * t1633 + t1544 * t1566;
t1648 = t1500 * t1508;
t1567 = cos(qJ(3,1));
t1632 = t1557 * t1564;
t1501 = -pkin(1) * t1632 + t1544 * t1567;
t1647 = t1501 * t1510;
t1502 = t1544 * t1562 + t1565 * t1670;
t1646 = t1502 * t1505;
t1645 = t1502 * t1506;
t1503 = t1544 * t1563 + t1566 * t1670;
t1644 = t1503 * t1507;
t1643 = t1503 * t1508;
t1504 = t1544 * t1564 + t1567 * t1670;
t1642 = t1504 * t1509;
t1641 = t1504 * t1510;
t1640 = t1505 / 0.2e1;
t1638 = t1507 / 0.2e1;
t1636 = t1509 / 0.2e1;
t1551 = sin(t1559);
t1552 = sin(t1560);
t1553 = sin(t1561);
t1554 = cos(t1559);
t1555 = cos(t1560);
t1556 = cos(t1561);
t1574 = t1490 * t1493 * t1506 / 0.4e1 + t1491 * t1494 * t1508 / 0.4e1 + t1492 * t1495 * t1510 / 0.4e1;
t1626 = (t1551 * t1554 + t1552 * t1555 + t1553 * t1556 + t1571 * t1574) * MDP(4) + MDP(1) * t1574;
t1625 = t1671 * (t1490 * t1673 + t1491 * t1674 + t1492 * t1675);
t1624 = t1671 * (t1493 * t1673 + t1494 * t1674 + t1495 * t1675);
t1478 = -pkin(2) * t1565 + (-t1558 * t1565 + t1634) * pkin(1);
t1619 = t1478 * t1655;
t1618 = t1478 * t1652;
t1617 = t1478 * t1496;
t1479 = -pkin(2) * t1566 + (-t1558 * t1566 + t1633) * pkin(1);
t1616 = t1479 * t1654;
t1615 = t1479 * t1651;
t1614 = t1479 * t1497;
t1480 = pkin(2) * t1567 + (t1558 * t1567 - t1632) * pkin(1);
t1613 = t1480 * t1653;
t1612 = t1480 * t1650;
t1611 = t1480 * t1498;
t1610 = t1490 * t1649;
t1609 = t1490 * t1646;
t1608 = t1490 * t1645;
t1607 = t1491 * t1648;
t1606 = t1491 * t1644;
t1605 = t1491 * t1643;
t1604 = t1492 * t1647;
t1603 = t1492 * t1642;
t1602 = t1492 * t1641;
t1601 = t1493 * t1649;
t1600 = t1493 * t1646;
t1599 = t1493 * t1645;
t1598 = t1494 * t1648;
t1597 = t1494 * t1644;
t1596 = t1494 * t1643;
t1595 = t1495 * t1647;
t1594 = t1495 * t1642;
t1593 = t1495 * t1641;
t1592 = t1499 * t1639;
t1591 = t1500 * t1637;
t1590 = t1501 * t1635;
t1589 = t1502 * t1496;
t1588 = t1502 * t1639;
t1587 = t1503 * t1497;
t1586 = t1503 * t1637;
t1585 = t1504 * t1498;
t1584 = t1504 * t1635;
t1583 = t1505 * t1631;
t1582 = t1507 * t1631;
t1581 = t1509 * t1631;
t1580 = t1454 * t1583;
t1579 = t1455 * t1582;
t1578 = t1456 * t1581;
t1474 = t1495 * t1636;
t1473 = t1494 * t1638;
t1472 = t1493 * t1640;
t1471 = t1492 * t1636;
t1470 = t1491 * t1638;
t1469 = t1490 * t1640;
t1467 = 0.2e1 * t1498 + t1477;
t1465 = 0.2e1 * t1497 + t1476;
t1463 = 0.2e1 * t1496 + t1475;
t1462 = t1498 + t1477 / 0.2e1;
t1461 = t1497 + t1476 / 0.2e1;
t1460 = t1496 + t1475 / 0.2e1;
t1453 = t1459 * t1581;
t1452 = t1458 * t1582;
t1451 = t1457 * t1583;
t1450 = t1474 + t1453;
t1449 = 0.2e1 * t1474 + t1453;
t1448 = t1473 + t1452;
t1447 = 0.2e1 * t1473 + t1452;
t1446 = t1472 + t1451;
t1445 = 0.2e1 * t1472 + t1451;
t1444 = t1471 - t1578;
t1443 = 0.2e1 * t1471 - t1578;
t1442 = t1470 - t1579;
t1441 = 0.2e1 * t1470 - t1579;
t1440 = t1469 - t1580;
t1439 = 0.2e1 * t1469 - t1580;
t1438 = t1474 + t1453 / 0.2e1;
t1437 = t1473 + t1452 / 0.2e1;
t1436 = t1472 + t1451 / 0.2e1;
t1435 = t1471 - t1578 / 0.2e1;
t1434 = t1470 - t1579 / 0.2e1;
t1433 = t1469 - t1580 / 0.2e1;
t1 = [(t1551 ^ 2 + t1552 ^ 2 + t1553 ^ 2) * MDP(4) + MDP(8) + ((-t1454 * t1610 - t1455 * t1607 - t1456 * t1604) * MDP(6) + (t1454 * t1608 + t1455 * t1605 + t1456 * t1602) * MDP(7)) * t1630 + (t1490 ^ 2 * t1506 + t1491 ^ 2 * t1508 + t1492 ^ 2 * t1510) * t1680 + (t1440 * t1655 + t1442 * t1654 + t1444 * t1653) * t1679 + (-t1439 * t1619 - t1441 * t1616 + t1443 * t1613) * t1678 + (-t1433 * t1609 - t1434 * t1606 - t1435 * t1603) * t1676 + (-t1440 * t1667 - t1442 * t1666 - t1444 * t1665) * t1672; ((-t1454 * t1601 - t1455 * t1598 - t1456 * t1595) * MDP(6) + (t1454 * t1599 + t1455 * t1596 + t1456 * t1593) * MDP(7)) * t1630 + (t1446 * t1655 + t1448 * t1654 + t1450 * t1653) * t1679 + (-t1445 * t1619 - t1447 * t1616 + t1449 * t1613) * t1678 + (-t1436 * t1609 - t1437 * t1606 - t1438 * t1603) * t1676 + (-t1446 * t1667 - t1448 * t1666 - t1450 * t1665) * t1672 + t1626; (t1464 * t1655 + t1466 * t1654 + t1468 * t1653) * t1679 + (-t1463 * t1619 - t1465 * t1616 + t1467 * t1613) * t1678 + (-t1460 * t1609 - t1461 * t1606 - t1462 * t1603) * t1676 + ((-t1454 * t1661 - t1455 * t1660 - t1456 * t1659) * MDP(5) + (-t1454 * t1592 - t1455 * t1591 - t1456 * t1590) * MDP(6) + (t1454 * t1588 + t1455 * t1586 + t1456 * t1584) * MDP(7)) * t1631 + t1625; ((t1457 * t1610 + t1458 * t1607 + t1459 * t1604) * MDP(6) + (-t1457 * t1608 - t1458 * t1605 - t1459 * t1602) * MDP(7)) * t1630 + (t1440 * t1652 + t1442 * t1651 + t1444 * t1650) * t1679 + (-t1439 * t1618 - t1441 * t1615 + t1443 * t1612) * t1678 + (-t1433 * t1600 - t1434 * t1597 - t1435 * t1594) * t1676 + (t1440 * t1664 + t1442 * t1663 + t1444 * t1662) * t1672 + t1626; (t1554 ^ 2 + t1555 ^ 2 + t1556 ^ 2) * MDP(4) + MDP(8) + ((t1457 * t1601 + t1458 * t1598 + t1459 * t1595) * MDP(6) + (-t1457 * t1599 - t1458 * t1596 - t1459 * t1593) * MDP(7)) * t1630 + (t1493 ^ 2 * t1506 + t1494 ^ 2 * t1508 + t1495 ^ 2 * t1510) * t1680 + (t1446 * t1652 + t1448 * t1651 + t1450 * t1650) * t1679 + (-t1445 * t1618 - t1447 * t1615 + t1449 * t1612) * t1678 + (-t1436 * t1600 - t1437 * t1597 - t1438 * t1594) * t1676 + (t1446 * t1664 + t1448 * t1663 + t1450 * t1662) * t1672; (t1464 * t1652 + t1466 * t1651 + t1468 * t1650) * t1679 + (-t1463 * t1618 - t1465 * t1615 + t1467 * t1612) * t1678 + (-t1460 * t1600 - t1461 * t1597 - t1462 * t1594) * t1676 + ((t1457 * t1661 + t1458 * t1660 + t1459 * t1659) * MDP(5) + (t1457 * t1592 + t1458 * t1591 + t1459 * t1590) * MDP(6) + (-t1457 * t1588 - t1458 * t1586 - t1459 * t1584) * MDP(7)) * t1631 + t1624; (t1440 * t1496 + t1442 * t1497 + t1444 * t1498) * MDP(5) + (-t1439 * t1617 - t1441 * t1614 + t1443 * t1611) * MDP(6) + (-t1433 * t1589 - t1434 * t1587 - t1435 * t1585) * t1620 + ((t1440 * t1658 + t1442 * t1657 + t1444 * t1656) * MDP(5) + (t1487 * t1610 + t1488 * t1607 + t1489 * t1604) * t1678 + (-t1487 * t1608 - t1488 * t1605 - t1489 * t1602) * t1677) * t1570 + t1625; (t1446 * t1496 + t1448 * t1497 + t1450 * t1498) * MDP(5) + (-t1445 * t1617 - t1447 * t1614 + t1449 * t1611) * MDP(6) + (-t1436 * t1589 - t1437 * t1587 - t1438 * t1585) * t1620 + ((t1446 * t1658 + t1448 * t1657 + t1450 * t1656) * MDP(5) + (t1487 * t1601 + t1488 * t1598 + t1489 * t1595) * t1678 + (-t1487 * t1599 - t1488 * t1596 - t1489 * t1593) * t1677) * t1570 + t1624; (t1464 * t1496 + t1466 * t1497 + t1468 * t1498) * MDP(5) + (-t1463 * t1617 - t1465 * t1614 + t1467 * t1611) * MDP(6) + (-t1460 * t1589 - t1461 * t1587 - t1462 * t1585) * t1620 + MDP(8) + ((t1464 * t1658 + t1466 * t1657 + t1468 * t1656) * MDP(5) + (t1487 * t1592 + t1488 * t1591 + t1489 * t1590) * MDP(6) + (-t1487 * t1588 - t1488 * t1586 - t1489 * t1584) * MDP(7)) * t1570 + t1671 * (t1506 * t1541 ^ 2 + t1508 * t1542 ^ 2 + t1510 * t1543 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
