% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPRR1G3A0
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
%   see P3RPRR1G3A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RPRR1G3A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3A0_inertia_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3RPRR1G3A0_inertia_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:27:09
% EndTime: 2020-03-09 21:27:12
% DurationCPUTime: 2.48s
% Computational Cost: add. (5739->294), mult. (3564->461), div. (669->7), fcn. (2937->62), ass. (0->219)
t1583 = 0.1e1 / pkin(3);
t1653 = t1583 / 0.2e1;
t1584 = pkin(1) ^ 2;
t1702 = MDP(4) * t1584 + MDP(1);
t1711 = t1702 / 0.4e1;
t1710 = MDP(5) / 0.2e1;
t1709 = MDP(6) / 0.2e1;
t1708 = MDP(7) / 0.2e1;
t1652 = t1583 / 0.4e1;
t1642 = 0.2e1 * MDP(7);
t1707 = t1642 / 0.2e1;
t1577 = sin(qJ(3,1));
t1645 = pkin(7) + qJ(3,1);
t1529 = pkin(1) * sin(t1645) + t1577 * pkin(2);
t1526 = 0.1e1 / t1529 ^ 2;
t1559 = sin(qJ(1,1) + t1645);
t1660 = t1526 * t1559;
t1706 = -t1660 / 0.2e1;
t1576 = sin(qJ(3,2));
t1644 = pkin(7) + qJ(3,2);
t1528 = pkin(1) * sin(t1644) + t1576 * pkin(2);
t1524 = 0.1e1 / t1528 ^ 2;
t1558 = sin(qJ(1,2) + t1644);
t1662 = t1524 * t1558;
t1705 = -t1662 / 0.2e1;
t1575 = sin(qJ(3,3));
t1643 = pkin(7) + qJ(3,3);
t1527 = pkin(1) * sin(t1643) + t1575 * pkin(2);
t1522 = 0.1e1 / t1527 ^ 2;
t1557 = sin(qJ(1,3) + t1643);
t1664 = t1522 * t1557;
t1704 = -t1664 / 0.2e1;
t1703 = MDP(5) * t1653;
t1572 = legFrame(3,2);
t1649 = qJ(1,3) + pkin(7);
t1554 = t1572 + t1649;
t1548 = qJ(3,3) + t1554;
t1600 = -t1572 + t1649;
t1599 = qJ(3,3) + t1600;
t1701 = sin(t1548) - sin(t1599);
t1573 = legFrame(2,2);
t1650 = qJ(1,2) + pkin(7);
t1555 = t1573 + t1650;
t1549 = qJ(3,2) + t1555;
t1601 = -t1573 + t1650;
t1598 = qJ(3,2) + t1601;
t1700 = sin(t1549) - sin(t1598);
t1574 = legFrame(1,2);
t1651 = qJ(1,1) + pkin(7);
t1556 = t1574 + t1651;
t1550 = qJ(3,1) + t1556;
t1602 = -t1574 + t1651;
t1597 = qJ(3,1) + t1602;
t1699 = sin(t1550) - sin(t1597);
t1570 = sin(pkin(7));
t1698 = pkin(1) * t1570;
t1695 = qJ(1,1) - t1574;
t1694 = qJ(1,2) - t1573;
t1693 = qJ(1,3) - t1572;
t1561 = qJ(1,3) + t1572;
t1473 = t1701 * pkin(3) + (sin(t1554) - sin(t1600)) * pkin(2) + (sin(t1561) - sin(t1693)) * pkin(1);
t1521 = 0.1e1 / t1527;
t1692 = t1473 * t1521;
t1562 = qJ(1,2) + t1573;
t1474 = t1700 * pkin(3) + (sin(t1555) - sin(t1601)) * pkin(2) + (sin(t1562) - sin(t1694)) * pkin(1);
t1523 = 0.1e1 / t1528;
t1691 = t1474 * t1523;
t1563 = qJ(1,1) + t1574;
t1475 = t1699 * pkin(3) + (sin(t1556) - sin(t1602)) * pkin(2) + (sin(t1563) - sin(t1695)) * pkin(1);
t1525 = 0.1e1 / t1529;
t1690 = t1475 * t1525;
t1512 = cos(t1599) + cos(t1548);
t1476 = -t1512 * pkin(3) + (-cos(t1600) - cos(t1554)) * pkin(2) + (-cos(t1693) - cos(t1561)) * pkin(1);
t1689 = t1476 * t1521;
t1513 = cos(t1598) + cos(t1549);
t1477 = -t1513 * pkin(3) + (-cos(t1601) - cos(t1555)) * pkin(2) + (-cos(t1694) - cos(t1562)) * pkin(1);
t1688 = t1477 * t1523;
t1514 = cos(t1597) + cos(t1550);
t1478 = -t1514 * pkin(3) + (-cos(t1602) - cos(t1556)) * pkin(2) + (-cos(t1695) - cos(t1563)) * pkin(1);
t1687 = t1478 * t1525;
t1506 = pkin(3) * t1557 + pkin(2) * sin(t1649) + sin(qJ(1,3)) * pkin(1);
t1683 = t1506 * t1521;
t1494 = t1583 * t1683;
t1659 = t1557 * t1521;
t1482 = t1494 - t1659;
t1686 = t1482 * t1521;
t1507 = pkin(3) * t1558 + pkin(2) * sin(t1650) + sin(qJ(1,2)) * pkin(1);
t1682 = t1507 * t1523;
t1495 = t1583 * t1682;
t1658 = t1558 * t1523;
t1484 = t1495 - t1658;
t1685 = t1484 * t1523;
t1508 = pkin(3) * t1559 + pkin(2) * sin(t1651) + sin(qJ(1,1)) * pkin(1);
t1681 = t1508 * t1525;
t1496 = t1583 * t1681;
t1657 = t1559 * t1525;
t1486 = t1496 - t1657;
t1684 = t1486 * t1525;
t1680 = t1701 * t1521;
t1679 = t1700 * t1523;
t1678 = t1699 * t1525;
t1677 = t1512 * t1521;
t1676 = t1513 * t1523;
t1675 = t1514 * t1525;
t1571 = cos(pkin(7));
t1560 = pkin(1) * t1571 + pkin(2);
t1578 = cos(qJ(3,3));
t1656 = t1570 * t1575;
t1515 = -pkin(1) * t1656 + t1560 * t1578;
t1674 = t1515 * t1522;
t1579 = cos(qJ(3,2));
t1655 = t1570 * t1576;
t1516 = -pkin(1) * t1655 + t1560 * t1579;
t1673 = t1516 * t1524;
t1580 = cos(qJ(3,1));
t1654 = t1570 * t1577;
t1517 = -pkin(1) * t1654 + t1560 * t1580;
t1672 = t1517 * t1526;
t1518 = t1560 * t1575 + t1578 * t1698;
t1671 = t1518 * t1521;
t1670 = t1518 * t1522;
t1519 = t1560 * t1576 + t1579 * t1698;
t1669 = t1519 * t1523;
t1668 = t1519 * t1524;
t1520 = t1560 * t1577 + t1580 * t1698;
t1667 = t1520 * t1525;
t1666 = t1520 * t1526;
t1665 = t1521 / 0.2e1;
t1663 = t1523 / 0.2e1;
t1661 = t1525 / 0.2e1;
t1564 = sin(t1572);
t1565 = sin(t1573);
t1566 = sin(t1574);
t1567 = cos(t1572);
t1568 = cos(t1573);
t1569 = cos(t1574);
t1587 = -t1701 * t1512 * t1522 / 0.4e1 - t1700 * t1513 * t1524 / 0.4e1 - t1699 * t1514 * t1526 / 0.4e1;
t1648 = MDP(1) * t1587 + (t1564 * t1567 + t1565 * t1568 + t1566 * t1569 + t1584 * t1587) * MDP(4);
t1647 = t1702 * (-t1699 * t1706 - t1700 * t1705 - t1701 * t1704);
t1646 = t1702 * (t1512 * t1704 + t1513 * t1705 + t1514 * t1706);
t1497 = -t1579 * pkin(2) + (-t1571 * t1579 + t1655) * pkin(1);
t1641 = t1497 * t1679;
t1640 = t1497 * t1676;
t1639 = t1497 * t1658;
t1498 = t1578 * pkin(2) + (t1571 * t1578 - t1656) * pkin(1);
t1638 = t1498 * t1680;
t1637 = t1498 * t1677;
t1636 = t1498 * t1659;
t1499 = t1580 * pkin(2) + (t1571 * t1580 - t1654) * pkin(1);
t1635 = t1499 * t1678;
t1634 = t1499 * t1675;
t1633 = t1499 * t1657;
t1632 = t1701 * t1674;
t1631 = t1701 * t1671;
t1630 = t1701 * t1670;
t1629 = t1700 * t1673;
t1628 = t1700 * t1669;
t1627 = t1700 * t1668;
t1626 = t1699 * t1672;
t1625 = t1699 * t1667;
t1624 = t1699 * t1666;
t1623 = t1512 * t1674;
t1622 = t1512 * t1671;
t1621 = t1512 * t1670;
t1620 = t1513 * t1673;
t1619 = t1513 * t1669;
t1618 = t1513 * t1668;
t1617 = t1514 * t1672;
t1616 = t1514 * t1667;
t1615 = t1514 * t1666;
t1614 = t1515 * t1664;
t1613 = t1516 * t1662;
t1612 = t1517 * t1660;
t1611 = t1518 * t1659;
t1610 = t1518 * t1664;
t1609 = t1519 * t1658;
t1608 = t1519 * t1662;
t1607 = t1520 * t1657;
t1606 = t1520 * t1660;
t1605 = t1521 * t1653;
t1604 = t1523 * t1653;
t1603 = t1525 * t1653;
t1596 = t1473 * t1605;
t1595 = t1474 * t1604;
t1594 = t1475 * t1603;
t1493 = t1514 * t1661;
t1492 = t1513 * t1663;
t1491 = t1512 * t1665;
t1490 = t1699 * t1661;
t1489 = t1700 * t1663;
t1488 = t1701 * t1665;
t1487 = t1496 - 0.2e1 * t1657;
t1485 = t1495 - 0.2e1 * t1658;
t1483 = t1494 - 0.2e1 * t1659;
t1481 = -t1657 + t1496 / 0.2e1;
t1480 = -t1658 + t1495 / 0.2e1;
t1479 = -t1659 + t1494 / 0.2e1;
t1472 = t1478 * t1603;
t1471 = t1477 * t1604;
t1470 = t1476 * t1605;
t1469 = t1493 + t1472;
t1468 = 0.2e1 * t1493 + t1472;
t1467 = t1492 + t1471;
t1466 = 0.2e1 * t1492 + t1471;
t1465 = t1491 + t1470;
t1464 = 0.2e1 * t1491 + t1470;
t1463 = -t1490 + t1594;
t1462 = -0.2e1 * t1490 + t1594;
t1461 = -t1489 + t1595;
t1460 = -0.2e1 * t1489 + t1595;
t1459 = -t1488 + t1596;
t1458 = -0.2e1 * t1488 + t1596;
t1457 = t1493 + t1472 / 0.2e1;
t1456 = t1492 + t1471 / 0.2e1;
t1455 = t1491 + t1470 / 0.2e1;
t1454 = -t1490 + t1594 / 0.2e1;
t1453 = -t1489 + t1595 / 0.2e1;
t1452 = -t1488 + t1596 / 0.2e1;
t1 = [(t1564 ^ 2 + t1565 ^ 2 + t1566 ^ 2) * MDP(4) + MDP(8) + ((t1476 * t1623 + t1477 * t1620 + t1478 * t1617) * MDP(6) + (-t1476 * t1621 - t1477 * t1618 - t1478 * t1615) * MDP(7)) * t1652 + (t1512 ^ 2 * t1522 + t1513 ^ 2 * t1524 + t1514 ^ 2 * t1526) * t1711 + (t1465 * t1677 + t1467 * t1676 + t1469 * t1675) * t1710 + (t1464 * t1637 - t1466 * t1640 + t1468 * t1634) * t1709 + (-t1455 * t1622 - t1456 * t1619 - t1457 * t1616) * t1707 + (t1465 * t1689 + t1467 * t1688 + t1469 * t1687) * t1703; ((-t1476 * t1632 - t1477 * t1629 - t1478 * t1626) * MDP(6) + (t1476 * t1630 + t1477 * t1627 + t1478 * t1624) * MDP(7)) * t1652 + (t1459 * t1677 + t1461 * t1676 + t1463 * t1675) * t1710 + (t1458 * t1637 - t1460 * t1640 + t1462 * t1634) * t1709 + (-t1452 * t1622 - t1453 * t1619 - t1454 * t1616) * t1707 + (t1459 * t1689 + t1461 * t1688 + t1463 * t1687) * t1703 + t1648; (t1482 * t1677 + t1484 * t1676 + t1486 * t1675) * t1710 + (t1483 * t1637 - t1485 * t1640 + t1487 * t1634) * t1709 + (-t1479 * t1622 - t1480 * t1619 - t1481 * t1616) * t1707 + ((t1476 * t1686 + t1477 * t1685 + t1478 * t1684) * MDP(5) + (-t1476 * t1614 - t1477 * t1613 - t1478 * t1612) * MDP(6) + (t1476 * t1610 + t1477 * t1608 + t1478 * t1606) * MDP(7)) * t1653 + t1646; ((t1473 * t1623 + t1474 * t1620 + t1475 * t1617) * MDP(6) + (-t1473 * t1621 - t1474 * t1618 - t1475 * t1615) * MDP(7)) * t1652 + (-t1465 * t1680 - t1467 * t1679 - t1469 * t1678) * t1710 + (-t1464 * t1638 + t1466 * t1641 - t1468 * t1635) * t1709 + (t1455 * t1631 + t1456 * t1628 + t1457 * t1625) * t1707 + (t1465 * t1692 + t1467 * t1691 + t1469 * t1690) * t1703 + t1648; (t1567 ^ 2 + t1568 ^ 2 + t1569 ^ 2) * MDP(4) + MDP(8) + ((-t1473 * t1632 - t1474 * t1629 - t1475 * t1626) * MDP(6) + (t1473 * t1630 + t1474 * t1627 + t1475 * t1624) * MDP(7)) * t1652 + (t1522 * t1701 ^ 2 + t1524 * t1700 ^ 2 + t1526 * t1699 ^ 2) * t1711 + (-t1459 * t1680 - t1461 * t1679 - t1463 * t1678) * t1710 + (-t1458 * t1638 + t1460 * t1641 - t1462 * t1635) * t1709 + (t1452 * t1631 + t1453 * t1628 + t1454 * t1625) * t1707 + (t1459 * t1692 + t1461 * t1691 + t1463 * t1690) * t1703; (-t1482 * t1680 - t1484 * t1679 - t1486 * t1678) * t1710 + (-t1483 * t1638 + t1485 * t1641 - t1487 * t1635) * t1709 + (t1479 * t1631 + t1480 * t1628 + t1481 * t1625) * t1707 + ((t1473 * t1686 + t1474 * t1685 + t1475 * t1684) * MDP(5) + (-t1473 * t1614 - t1474 * t1613 - t1475 * t1612) * MDP(6) + (t1473 * t1610 + t1474 * t1608 + t1475 * t1606) * MDP(7)) * t1653 + t1647; (-t1465 * t1659 - t1467 * t1658 - t1469 * t1657) * MDP(5) + (-t1464 * t1636 + t1466 * t1639 - t1468 * t1633) * MDP(6) + (t1455 * t1611 + t1456 * t1609 + t1457 * t1607) * t1642 + ((t1465 * t1683 + t1467 * t1682 + t1469 * t1681) * MDP(5) + (t1506 * t1623 + t1507 * t1620 + t1508 * t1617) * t1709 + (-t1506 * t1621 - t1507 * t1618 - t1508 * t1615) * t1708) * t1583 + t1646; (-t1459 * t1659 - t1461 * t1658 - t1463 * t1657) * MDP(5) + (-t1458 * t1636 + t1460 * t1639 - t1462 * t1633) * MDP(6) + (t1452 * t1611 + t1453 * t1609 + t1454 * t1607) * t1642 + ((t1459 * t1683 + t1461 * t1682 + t1463 * t1681) * MDP(5) + (-t1506 * t1632 - t1507 * t1629 - t1508 * t1626) * t1709 + (t1506 * t1630 + t1507 * t1627 + t1508 * t1624) * t1708) * t1583 + t1647; (-t1482 * t1659 - t1484 * t1658 - t1486 * t1657) * MDP(5) + (-t1483 * t1636 + t1485 * t1639 - t1487 * t1633) * MDP(6) + (t1479 * t1611 + t1480 * t1609 + t1481 * t1607) * t1642 + MDP(8) + ((t1482 * t1683 + t1484 * t1682 + t1486 * t1681) * MDP(5) + (-t1506 * t1614 - t1507 * t1613 - t1508 * t1612) * MDP(6) + (t1506 * t1610 + t1507 * t1608 + t1508 * t1606) * MDP(7)) * t1583 + t1702 * (t1522 * t1557 ^ 2 + t1524 * t1558 ^ 2 + t1526 * t1559 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
