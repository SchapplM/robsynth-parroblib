% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [10x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRR1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRR1G1P1A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(10,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1G1P1A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRR1G1P1A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1G1P1A0_coriolisvec_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1G1P1A0_coriolisvec_para_pf_mdp: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1G1P1A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1G1P1A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'P3RRR1G1P1A0_coriolisvec_para_pf_mdp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:43
% EndTime: 2019-05-03 15:38:47
% DurationCPUTime: 4.65s
% Computational Cost: add. (24382->253), mult. (38835->510), div. (2856->9), fcn. (29350->26), ass. (0->213)
t1558 = legFrame(3,3);
t1546 = sin(t1558);
t1549 = cos(t1558);
t1562 = sin(qJ(1,3));
t1568 = cos(qJ(1,3));
t1573 = xDP(3);
t1574 = xDP(2);
t1575 = xDP(1);
t1577 = koppelP(3,2);
t1580 = koppelP(3,1);
t1528 = -t1580 * t1562 + t1568 * t1577;
t1531 = t1562 * t1577 + t1568 * t1580;
t1576 = xP(3);
t1552 = sin(t1576);
t1553 = cos(t1576);
t1691 = (t1528 * t1553 + t1531 * t1552) * t1549 - t1546 * (-t1528 * t1552 + t1531 * t1553);
t1522 = t1552 * t1580 + t1553 * t1577;
t1495 = t1522 * t1573 - t1575;
t1525 = -t1552 * t1577 + t1553 * t1580;
t1498 = t1525 * t1573 + t1574;
t1555 = qJ(1,3) + qJ(2,3);
t1540 = sin(t1555);
t1543 = cos(t1555);
t1694 = t1543 * (t1495 * t1549 - t1498 * t1546) - (t1495 * t1546 + t1498 * t1549) * t1540;
t1462 = pkin(1) * ((-t1562 * t1574 - t1568 * t1575) * t1549 - t1546 * (-t1562 * t1575 + t1568 * t1574) + t1691 * t1573) + t1694 * pkin(2);
t1519 = t1540 * t1568 - t1543 * t1562;
t1688 = 0.1e1 / t1519;
t1678 = t1462 * t1688;
t1559 = legFrame(2,3);
t1547 = sin(t1559);
t1550 = cos(t1559);
t1564 = sin(qJ(1,2));
t1570 = cos(qJ(1,2));
t1578 = koppelP(2,2);
t1581 = koppelP(2,1);
t1529 = -t1581 * t1564 + t1570 * t1578;
t1532 = t1564 * t1578 + t1570 * t1581;
t1690 = (t1529 * t1553 + t1532 * t1552) * t1550 - t1547 * (-t1529 * t1552 + t1532 * t1553);
t1523 = t1552 * t1581 + t1553 * t1578;
t1496 = t1523 * t1573 - t1575;
t1526 = -t1552 * t1578 + t1553 * t1581;
t1499 = t1526 * t1573 + t1574;
t1556 = qJ(1,2) + qJ(2,2);
t1541 = sin(t1556);
t1544 = cos(t1556);
t1693 = t1544 * (t1496 * t1550 - t1499 * t1547) - (t1496 * t1547 + t1499 * t1550) * t1541;
t1463 = pkin(1) * ((-t1564 * t1574 - t1570 * t1575) * t1550 - t1547 * (-t1564 * t1575 + t1570 * t1574) + t1690 * t1573) + t1693 * pkin(2);
t1520 = t1541 * t1570 - t1544 * t1564;
t1687 = 0.1e1 / t1520;
t1677 = t1463 * t1687;
t1560 = legFrame(1,3);
t1548 = sin(t1560);
t1551 = cos(t1560);
t1566 = sin(qJ(1,1));
t1572 = cos(qJ(1,1));
t1579 = koppelP(1,2);
t1582 = koppelP(1,1);
t1530 = -t1582 * t1566 + t1572 * t1579;
t1533 = t1566 * t1579 + t1572 * t1582;
t1689 = (t1530 * t1553 + t1533 * t1552) * t1551 - t1548 * (-t1530 * t1552 + t1533 * t1553);
t1524 = t1552 * t1582 + t1553 * t1579;
t1497 = t1524 * t1573 - t1575;
t1527 = -t1552 * t1579 + t1553 * t1582;
t1500 = t1527 * t1573 + t1574;
t1557 = qJ(1,1) + qJ(2,1);
t1542 = sin(t1557);
t1545 = cos(t1557);
t1692 = t1545 * (t1497 * t1551 - t1500 * t1548) - (t1497 * t1548 + t1500 * t1551) * t1542;
t1464 = pkin(1) * ((-t1566 * t1574 - t1572 * t1575) * t1551 - t1548 * (-t1566 * t1575 + t1572 * t1574) + t1689 * t1573) + t1692 * pkin(2);
t1521 = t1542 * t1572 - t1545 * t1566;
t1686 = 0.1e1 / t1521;
t1676 = t1464 * t1686;
t1685 = 0.2e1 * pkin(2);
t1584 = 0.1e1 / pkin(2);
t1585 = 0.1e1 / pkin(1);
t1653 = t1584 * t1585;
t1459 = t1653 * t1678;
t1675 = t1694 * t1688;
t1468 = t1585 * t1675;
t1454 = -t1468 + t1459;
t1598 = t1540 * t1562 + t1543 * t1568;
t1447 = -pkin(2) * t1454 + t1598 * t1675;
t1511 = 0.1e1 / t1519 ^ 2;
t1586 = 0.1e1 / pkin(1) ^ 2;
t1501 = t1540 * t1549 + t1543 * t1546;
t1483 = pkin(1) * (t1546 * t1568 + t1549 * t1562) + t1501 * pkin(2);
t1502 = -t1540 * t1546 + t1543 * t1549;
t1486 = pkin(1) * (-t1546 * t1562 + t1549 * t1568) + t1502 * pkin(2);
t1589 = (-t1483 * t1522 - t1486 * t1525) * t1584 * t1688;
t1604 = -t1501 * t1522 - t1502 * t1525;
t1592 = t1604 * t1688;
t1637 = t1688 * t1678;
t1628 = t1454 * t1637;
t1629 = t1454 * t1462 * (pkin(1) * t1598 + pkin(2)) * t1584;
t1554 = t1573 ^ 2;
t1654 = t1554 * t1585;
t1583 = pkin(2) ^ 2;
t1681 = (t1454 * t1583 + (-t1675 + t1598 * (-t1468 + t1459 / 0.2e1) * t1685) * pkin(1)) * t1584;
t1438 = (t1628 + (-t1629 - (-t1447 - t1681) * t1694) * t1511) * t1586 + (-t1589 + t1592) * t1654;
t1684 = t1438 * t1688;
t1460 = t1653 * t1677;
t1674 = t1693 * t1687;
t1469 = t1585 * t1674;
t1456 = -t1469 + t1460;
t1597 = t1541 * t1564 + t1544 * t1570;
t1448 = -pkin(2) * t1456 + t1597 * t1674;
t1513 = 0.1e1 / t1520 ^ 2;
t1503 = t1541 * t1550 + t1544 * t1547;
t1484 = pkin(1) * (t1547 * t1570 + t1550 * t1564) + t1503 * pkin(2);
t1504 = -t1541 * t1547 + t1544 * t1550;
t1487 = pkin(1) * (-t1547 * t1564 + t1550 * t1570) + t1504 * pkin(2);
t1588 = (-t1484 * t1523 - t1487 * t1526) * t1584 * t1687;
t1603 = -t1503 * t1523 - t1504 * t1526;
t1591 = t1603 * t1687;
t1636 = t1687 * t1677;
t1624 = t1456 * t1636;
t1625 = t1456 * t1463 * (pkin(1) * t1597 + pkin(2)) * t1584;
t1680 = (t1456 * t1583 + (-t1674 + t1597 * (-t1469 + t1460 / 0.2e1) * t1685) * pkin(1)) * t1584;
t1439 = (t1624 + (-t1625 - (-t1448 - t1680) * t1693) * t1513) * t1586 + (-t1588 + t1591) * t1654;
t1683 = t1439 * t1687;
t1461 = t1653 * t1676;
t1673 = t1692 * t1686;
t1470 = t1585 * t1673;
t1458 = -t1470 + t1461;
t1596 = t1542 * t1566 + t1545 * t1572;
t1449 = -pkin(2) * t1458 + t1596 * t1673;
t1515 = 0.1e1 / t1521 ^ 2;
t1505 = t1542 * t1551 + t1545 * t1548;
t1485 = pkin(1) * (t1548 * t1572 + t1551 * t1566) + t1505 * pkin(2);
t1506 = -t1542 * t1548 + t1545 * t1551;
t1488 = pkin(1) * (-t1548 * t1566 + t1551 * t1572) + t1506 * pkin(2);
t1587 = (-t1485 * t1524 - t1488 * t1527) * t1584 * t1686;
t1602 = -t1505 * t1524 - t1506 * t1527;
t1590 = t1602 * t1686;
t1635 = t1686 * t1676;
t1620 = t1458 * t1635;
t1621 = t1458 * t1464 * (pkin(1) * t1596 + pkin(2)) * t1584;
t1679 = (t1458 * t1583 + (-t1673 + t1596 * (-t1470 + t1461 / 0.2e1) * t1685) * pkin(1)) * t1584;
t1440 = (t1620 + (-t1621 - (-t1449 - t1679) * t1692) * t1515) * t1586 + (-t1587 + t1590) * t1654;
t1682 = t1440 * t1686;
t1477 = (t1522 * t1549 - t1525 * t1546) * t1543 - (t1522 * t1546 + t1525 * t1549) * t1540;
t1672 = t1477 * t1688;
t1478 = (t1523 * t1550 - t1526 * t1547) * t1544 - (t1523 * t1547 + t1526 * t1550) * t1541;
t1671 = t1478 * t1687;
t1479 = (t1524 * t1551 - t1527 * t1548) * t1545 - (t1524 * t1548 + t1527 * t1551) * t1542;
t1670 = t1479 * t1686;
t1666 = t1501 * t1688;
t1665 = t1502 * t1688;
t1664 = t1503 * t1687;
t1663 = t1504 * t1687;
t1662 = t1505 * t1686;
t1661 = t1506 * t1686;
t1561 = sin(qJ(2,3));
t1660 = t1688 * t1561;
t1567 = cos(qJ(2,3));
t1659 = t1688 * t1567;
t1563 = sin(qJ(2,2));
t1658 = t1687 * t1563;
t1569 = cos(qJ(2,2));
t1657 = t1687 * t1569;
t1565 = sin(qJ(2,1));
t1656 = t1686 * t1565;
t1571 = cos(qJ(2,1));
t1655 = t1686 * t1571;
t1435 = (0.2e1 * t1628 + (-t1629 - (-0.2e1 * t1447 - t1681) * t1694) * t1511) * t1586 + (-t1589 + 0.2e1 * t1592) * t1654;
t1652 = t1435 * t1666;
t1651 = t1435 * t1660;
t1650 = t1435 * t1659;
t1436 = (0.2e1 * t1624 + (-t1625 - (-0.2e1 * t1448 - t1680) * t1693) * t1513) * t1586 + (-t1588 + 0.2e1 * t1591) * t1654;
t1649 = t1436 * t1664;
t1648 = t1436 * t1658;
t1647 = t1436 * t1657;
t1437 = (0.2e1 * t1620 + (-t1621 - (-0.2e1 * t1449 - t1679) * t1692) * t1515) * t1586 + (-t1587 + 0.2e1 * t1590) * t1654;
t1646 = t1437 * t1662;
t1645 = t1437 * t1656;
t1644 = t1437 * t1655;
t1441 = t1447 * t1694 * t1511 * t1586 + (t1454 * t1586 * t1678 + t1604 * t1654) * t1688;
t1643 = t1441 * t1660;
t1642 = t1441 * t1659;
t1442 = t1448 * t1693 * t1513 * t1586 + (t1456 * t1586 * t1677 + t1603 * t1654) * t1687;
t1641 = t1442 * t1658;
t1640 = t1442 * t1657;
t1443 = t1449 * t1692 * t1515 * t1586 + (t1458 * t1586 * t1676 + t1602 * t1654) * t1686;
t1639 = t1443 * t1656;
t1638 = t1443 * t1655;
t1634 = t1694 ^ 2 * t1511 * t1688;
t1633 = t1693 ^ 2 * t1513 * t1687;
t1632 = t1692 ^ 2 * t1515 * t1686;
t1453 = -0.2e1 * t1468 + t1459;
t1631 = t1453 * t1462 * t1501 * t1511;
t1630 = t1453 * t1637;
t1455 = -0.2e1 * t1469 + t1460;
t1627 = t1455 * t1463 * t1503 * t1513;
t1626 = t1455 * t1636;
t1457 = -0.2e1 * t1470 + t1461;
t1623 = t1457 * t1464 * t1505 * t1515;
t1622 = t1457 * t1635;
t1619 = t1561 * t1634;
t1618 = t1567 * t1634;
t1617 = t1563 * t1633;
t1616 = t1569 * t1633;
t1615 = t1565 * t1632;
t1614 = t1571 * t1632;
t1613 = t1561 * t1630;
t1612 = t1567 * t1630;
t1611 = t1563 * t1626;
t1610 = t1569 * t1626;
t1609 = t1565 * t1622;
t1608 = t1571 * t1622;
t1467 = -pkin(1) * t1689 - t1479 * pkin(2);
t1466 = -pkin(1) * t1690 - t1478 * pkin(2);
t1465 = -pkin(1) * t1691 - t1477 * pkin(2);
t1 = [(t1502 * t1650 + t1504 * t1647 + t1506 * t1644) * MDP(5) + (-t1502 * t1651 - t1504 * t1648 - t1506 * t1645) * MDP(6) + (-MDP(8) * t1553 + MDP(9) * t1552) * t1554 + ((t1441 * t1665 + t1442 * t1663 + t1443 * t1661) * MDP(1) + (t1438 * t1665 + t1439 * t1663 + t1440 * t1661) * MDP(4)) * t1585 + ((-t1486 * t1642 - t1487 * t1640 - t1488 * t1638) * MDP(5) + (t1486 * t1643 + t1487 * t1641 + t1488 * t1639) * MDP(6) + ((-t1486 * t1619 - t1487 * t1617 - t1488 * t1615) * MDP(5) + (-t1486 * t1618 - t1487 * t1616 - t1488 * t1614) * MDP(6)) * t1586 + ((-t1486 * t1684 - t1487 * t1683 - t1488 * t1682) * MDP(4) + (-t1502 * t1613 - t1504 * t1611 - t1506 * t1609) * MDP(5) + (-t1502 * t1612 - t1504 * t1610 - t1506 * t1608) * MDP(6)) * t1585) * t1584; (t1567 * t1652 + t1569 * t1649 + t1571 * t1646) * MDP(5) + (-t1561 * t1652 - t1563 * t1649 - t1565 * t1646) * MDP(6) + (-MDP(8) * t1552 - MDP(9) * t1553) * t1554 + ((t1441 * t1666 + t1442 * t1664 + t1443 * t1662) * MDP(1) + (t1438 * t1666 + t1439 * t1664 + t1440 * t1662) * MDP(4)) * t1585 + ((-t1483 * t1642 - t1484 * t1640 - t1485 * t1638) * MDP(5) + (t1483 * t1643 + t1484 * t1641 + t1485 * t1639) * MDP(6) + ((-t1483 * t1619 - t1484 * t1617 - t1485 * t1615) * MDP(5) + (-t1483 * t1618 - t1484 * t1616 - t1485 * t1614) * MDP(6)) * t1586 + ((-t1483 * t1684 - t1484 * t1683 - t1485 * t1682) * MDP(4) + (-t1561 * t1631 - t1563 * t1627 - t1565 * t1623) * MDP(5) + (-t1567 * t1631 - t1569 * t1627 - t1571 * t1623) * MDP(6)) * t1585) * t1584; (-t1477 * t1650 - t1478 * t1647 - t1479 * t1644) * MDP(5) + (t1477 * t1651 + t1478 * t1648 + t1479 * t1645) * MDP(6) + ((-t1441 * t1672 - t1442 * t1671 - t1443 * t1670) * MDP(1) + (-t1438 * t1672 - t1439 * t1671 - t1440 * t1670) * MDP(4)) * t1585 + ((-t1465 * t1642 - t1466 * t1640 - t1467 * t1638) * MDP(5) + (t1465 * t1643 + t1466 * t1641 + t1467 * t1639) * MDP(6) + ((-t1465 * t1619 - t1466 * t1617 - t1467 * t1615) * MDP(5) + (-t1465 * t1618 - t1466 * t1616 - t1467 * t1614) * MDP(6)) * t1586 + ((-t1465 * t1684 - t1466 * t1683 - t1467 * t1682) * MDP(4) + (t1477 * t1613 + t1478 * t1611 + t1479 * t1609) * MDP(5) + (t1477 * t1612 + t1478 * t1610 + t1479 * t1608) * MDP(6)) * t1585) * t1584;];
taucX  = t1;
