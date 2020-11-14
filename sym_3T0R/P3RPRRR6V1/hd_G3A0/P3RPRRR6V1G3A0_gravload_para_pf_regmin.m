% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR6V1G3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR6V1G3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:42:32
% EndTime: 2020-08-06 18:42:33
% DurationCPUTime: 1.37s
% Computational Cost: add. (900->187), mult. (1068->291), div. (111->10), fcn. (957->53), ass. (0->163)
t1566 = sin(pkin(7));
t1584 = -pkin(6) - pkin(5);
t1526 = t1584 * t1566 - pkin(1);
t1578 = cos(qJ(1,3));
t1517 = t1526 * t1578;
t1567 = cos(pkin(7));
t1572 = sin(qJ(1,3));
t1630 = t1572 * t1584;
t1631 = t1572 * t1566;
t1677 = t1517 + pkin(2) * t1631 - (pkin(2) * t1578 - t1630) * t1567;
t1580 = cos(qJ(1,2));
t1518 = t1526 * t1580;
t1574 = sin(qJ(1,2));
t1626 = t1574 * t1584;
t1627 = t1574 * t1566;
t1676 = t1518 + pkin(2) * t1627 - (pkin(2) * t1580 - t1626) * t1567;
t1582 = cos(qJ(1,1));
t1519 = t1526 * t1582;
t1576 = sin(qJ(1,1));
t1622 = t1576 * t1584;
t1623 = t1576 * t1566;
t1675 = t1519 + pkin(2) * t1623 - (pkin(2) * t1582 - t1622) * t1567;
t1568 = legFrame(3,2);
t1545 = sin(t1568);
t1548 = cos(t1568);
t1514 = t1548 * g(1) - t1545 * g(2);
t1497 = -g(3) * t1572 + t1514 * t1578;
t1674 = t1497 / 0.2e1;
t1569 = legFrame(2,2);
t1546 = sin(t1569);
t1549 = cos(t1569);
t1515 = t1549 * g(1) - t1546 * g(2);
t1499 = -g(3) * t1574 + t1515 * t1580;
t1673 = t1499 / 0.2e1;
t1570 = legFrame(1,2);
t1547 = sin(t1570);
t1550 = cos(t1570);
t1516 = t1550 * g(1) - t1547 * g(2);
t1501 = -g(3) * t1576 + t1516 * t1582;
t1672 = t1501 / 0.2e1;
t1581 = cos(qJ(3,1));
t1537 = t1581 * pkin(3) + pkin(2);
t1544 = t1567 * pkin(1);
t1522 = 0.1e1 / (t1544 + t1537);
t1559 = qJ(1,1) + pkin(7);
t1540 = sin(t1559);
t1543 = cos(t1559);
t1649 = (g(3) * t1543 + t1516 * t1540) * t1522;
t1604 = t1581 * t1649;
t1671 = t1604 / 0.2e1;
t1575 = sin(qJ(3,1));
t1605 = t1575 * t1649;
t1670 = -t1605 / 0.2e1;
t1579 = cos(qJ(3,2));
t1536 = t1579 * pkin(3) + pkin(2);
t1521 = 0.1e1 / (t1544 + t1536);
t1558 = qJ(1,2) + pkin(7);
t1539 = sin(t1558);
t1542 = cos(t1558);
t1650 = (g(3) * t1542 + t1515 * t1539) * t1521;
t1606 = t1579 * t1650;
t1669 = t1606 / 0.2e1;
t1573 = sin(qJ(3,2));
t1607 = t1573 * t1650;
t1668 = -t1607 / 0.2e1;
t1577 = cos(qJ(3,3));
t1535 = t1577 * pkin(3) + pkin(2);
t1520 = 0.1e1 / (t1544 + t1535);
t1557 = qJ(1,3) + pkin(7);
t1538 = sin(t1557);
t1541 = cos(t1557);
t1651 = (g(3) * t1541 + t1514 * t1538) * t1520;
t1608 = t1577 * t1651;
t1667 = t1608 / 0.2e1;
t1571 = sin(qJ(3,3));
t1609 = t1571 * t1651;
t1666 = -t1609 / 0.2e1;
t1598 = g(3) * t1582 + t1516 * t1576;
t1532 = t1570 + t1559;
t1533 = -t1570 + t1559;
t1507 = cos(t1533) - cos(t1532);
t1643 = t1507 * t1522;
t1665 = t1598 * t1643;
t1504 = -sin(t1532) - sin(t1533);
t1646 = t1504 * t1522;
t1664 = t1598 * t1646;
t1599 = g(3) * t1580 + t1515 * t1574;
t1530 = t1569 + t1558;
t1531 = -t1569 + t1558;
t1506 = cos(t1531) - cos(t1530);
t1644 = t1506 * t1521;
t1663 = t1599 * t1644;
t1503 = -sin(t1530) - sin(t1531);
t1647 = t1503 * t1521;
t1662 = t1599 * t1647;
t1600 = g(3) * t1578 + t1514 * t1572;
t1528 = t1568 + t1557;
t1529 = -t1568 + t1557;
t1505 = cos(t1529) - cos(t1528);
t1645 = t1505 * t1520;
t1661 = t1600 * t1645;
t1502 = -sin(t1528) - sin(t1529);
t1648 = t1502 * t1520;
t1660 = t1600 * t1648;
t1638 = t1522 * t1543;
t1640 = t1521 * t1542;
t1642 = t1520 * t1541;
t1659 = -t1598 * t1638 - t1599 * t1640 - t1600 * t1642;
t1658 = 0.2e1 * pkin(1);
t1657 = 0.2e1 * pkin(2);
t1656 = 0.2e1 * t1584;
t1655 = pkin(1) / 0.2e1;
t1616 = pkin(7) + qJ(3,3);
t1619 = -pkin(7) + qJ(3,3);
t1654 = (t1541 * t1656 + t1572 * t1658 + t1538 * t1657 + (sin(qJ(1,3) - t1619) + sin(qJ(1,3) + t1616)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,3)) + t1571 * t1657 + (sin(t1616) + sin(t1619)) * pkin(1));
t1617 = pkin(7) + qJ(3,2);
t1620 = -pkin(7) + qJ(3,2);
t1653 = (t1542 * t1656 + t1574 * t1658 + t1539 * t1657 + (sin(qJ(1,2) - t1620) + sin(qJ(1,2) + t1617)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,2)) + t1573 * t1657 + (sin(t1617) + sin(t1620)) * pkin(1));
t1618 = pkin(7) + qJ(3,1);
t1621 = -pkin(7) + qJ(3,1);
t1652 = (t1543 * t1656 + t1576 * t1658 + t1540 * t1657 + (sin(qJ(1,1) - t1621) + sin(qJ(1,1) + t1618)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,1)) + t1575 * t1657 + (sin(t1618) + sin(t1621)) * pkin(1));
t1641 = t1520 / t1571;
t1639 = t1521 / t1573;
t1637 = t1522 / t1575;
t1636 = t1535 * t1578;
t1635 = t1536 * t1580;
t1634 = t1537 * t1582;
t1633 = t1571 * t1545;
t1632 = t1571 * t1548;
t1629 = t1573 * t1546;
t1628 = t1573 * t1549;
t1625 = t1575 * t1547;
t1624 = t1575 * t1550;
t1615 = pkin(3) * (-t1578 * t1567 + t1631) * t1577 ^ 2;
t1614 = pkin(3) * (-t1580 * t1567 + t1627) * t1579 ^ 2;
t1613 = pkin(3) * (-t1582 * t1567 + t1623) * t1581 ^ 2;
t1612 = ((-t1630 + t1636) * t1567 - t1517 - t1535 * t1631) * t1641;
t1611 = ((-t1626 + t1635) * t1567 - t1518 - t1536 * t1627) * t1639;
t1610 = ((-t1622 + t1634) * t1567 - t1519 - t1537 * t1623) * t1637;
t1511 = t1545 * g(1) + t1548 * g(2);
t1603 = t1511 * t1641;
t1512 = t1546 * g(1) + t1549 * g(2);
t1602 = t1512 * t1639;
t1513 = t1547 * g(1) + t1550 * g(2);
t1601 = t1513 * t1637;
t1597 = t1545 * t1612;
t1596 = t1548 * t1612;
t1595 = t1546 * t1611;
t1594 = t1549 * t1611;
t1593 = t1547 * t1610;
t1592 = t1550 * t1610;
t1591 = -g(3) * t1538 + t1514 * t1541;
t1590 = -g(3) * t1539 + t1515 * t1542;
t1589 = -g(3) * t1540 + t1516 * t1543;
t1585 = 0.1e1 / pkin(3);
t1534 = t1544 + pkin(2);
t1480 = t1513 * t1575 + t1589 * t1581;
t1479 = -t1513 * t1581 + t1589 * t1575;
t1478 = t1512 * t1573 + t1590 * t1579;
t1477 = -t1512 * t1579 + t1590 * t1573;
t1476 = t1511 * t1571 + t1591 * t1577;
t1475 = -t1511 * t1577 + t1591 * t1571;
t1 = [0, t1664 / 0.2e1 + t1662 / 0.2e1 + t1660 / 0.2e1, t1646 * t1672 + t1647 * t1673 + t1648 * t1674, -(-t1550 * t1613 + (pkin(3) * t1625 - t1675 * t1550) * t1581 + t1534 * t1625) * t1601 - (-t1549 * t1614 + (pkin(3) * t1629 - t1676 * t1549) * t1579 + t1534 * t1629) * t1602 - (-t1548 * t1615 + (pkin(3) * t1633 - t1677 * t1548) * t1577 + t1534 * t1633) * t1603 + (t1660 + t1662 + t1664) * t1655, 0, 0, 0, 0, 0, (-t1475 * t1596 - t1477 * t1594 - t1479 * t1592) * t1585 + t1502 * t1667 + t1503 * t1669 + t1504 * t1671, (-t1476 * t1596 - t1478 * t1594 - t1480 * t1592) * t1585 + t1502 * t1666 + t1503 * t1668 + t1504 * t1670, -g(1); 0, t1665 / 0.2e1 + t1663 / 0.2e1 + t1661 / 0.2e1, t1643 * t1672 + t1644 * t1673 + t1645 * t1674, -(t1547 * t1613 + (pkin(3) * t1624 + t1675 * t1547) * t1581 + t1534 * t1624) * t1601 - (t1546 * t1614 + (pkin(3) * t1628 + t1676 * t1546) * t1579 + t1534 * t1628) * t1602 - (t1545 * t1615 + (pkin(3) * t1632 + t1677 * t1545) * t1577 + t1534 * t1632) * t1603 + (t1661 + t1663 + t1665) * t1655, 0, 0, 0, 0, 0, (t1475 * t1597 + t1477 * t1595 + t1479 * t1593) * t1585 + t1505 * t1667 + t1506 * t1669 + t1507 * t1671, (t1476 * t1597 + t1478 * t1595 + t1480 * t1593) * t1585 + t1505 * t1666 + t1506 * t1668 + t1507 * t1670, -g(2); 0, t1659, -t1497 * t1642 - t1499 * t1640 - t1501 * t1638, t1581 * ((t1537 * t1576 + t1582 * t1584) * t1567 - t1526 * t1576 + t1566 * t1634) * t1601 + t1579 * ((t1536 * t1574 + t1580 * t1584) * t1567 - t1526 * t1574 + t1566 * t1635) * t1602 + t1577 * ((t1535 * t1572 + t1578 * t1584) * t1567 - t1526 * t1572 + t1566 * t1636) * t1603 + t1659 * pkin(1), 0, 0, 0, 0, 0, -t1541 * t1608 - t1542 * t1606 - t1543 * t1604 + (t1475 * t1654 + t1477 * t1653 + t1479 * t1652) * t1585, t1541 * t1609 + t1542 * t1607 + t1543 * t1605 + (t1476 * t1654 + t1478 * t1653 + t1480 * t1652) * t1585, -g(3);];
tau_reg  = t1;
