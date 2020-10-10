% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR12V1G3A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x15]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:11
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR12V1G3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:10:48
% EndTime: 2020-08-06 19:10:50
% DurationCPUTime: 1.53s
% Computational Cost: add. (1035->148), mult. (1884->282), div. (171->6), fcn. (1764->18), ass. (0->152)
t1647 = sin(qJ(2,3));
t1654 = cos(qJ(1,3));
t1707 = t1647 * t1654;
t1648 = sin(qJ(1,3));
t1755 = t1648 * pkin(4);
t1602 = qJ(3,3) * t1707 - t1755;
t1644 = legFrame(3,2);
t1629 = sin(t1644);
t1632 = cos(t1644);
t1653 = cos(qJ(2,3));
t1641 = t1653 ^ 2;
t1659 = pkin(1) + pkin(2);
t1701 = t1659 * t1647;
t1704 = t1654 * t1659;
t1746 = t1632 * qJ(3,3);
t1566 = (-t1629 * t1704 - t1746) * t1641 + (-t1629 * t1602 + t1632 * t1701) * t1653 + t1746;
t1649 = sin(qJ(2,2));
t1656 = cos(qJ(1,2));
t1706 = t1649 * t1656;
t1650 = sin(qJ(1,2));
t1754 = t1650 * pkin(4);
t1603 = qJ(3,2) * t1706 - t1754;
t1645 = legFrame(2,2);
t1630 = sin(t1645);
t1633 = cos(t1645);
t1655 = cos(qJ(2,2));
t1642 = t1655 ^ 2;
t1700 = t1659 * t1649;
t1703 = t1656 * t1659;
t1745 = t1633 * qJ(3,2);
t1567 = (-t1630 * t1703 - t1745) * t1642 + (-t1630 * t1603 + t1633 * t1700) * t1655 + t1745;
t1651 = sin(qJ(2,1));
t1658 = cos(qJ(1,1));
t1705 = t1651 * t1658;
t1652 = sin(qJ(1,1));
t1753 = t1652 * pkin(4);
t1604 = qJ(3,1) * t1705 - t1753;
t1646 = legFrame(1,2);
t1631 = sin(t1646);
t1634 = cos(t1646);
t1657 = cos(qJ(2,1));
t1643 = t1657 ^ 2;
t1699 = t1659 * t1651;
t1702 = t1658 * t1659;
t1744 = t1634 * qJ(3,1);
t1568 = (-t1631 * t1702 - t1744) * t1643 + (-t1631 * t1604 + t1634 * t1699) * t1657 + t1744;
t1626 = t1647 * qJ(3,3);
t1611 = t1659 * t1653 + t1626;
t1605 = 0.1e1 / t1611;
t1627 = t1649 * qJ(3,2);
t1612 = t1659 * t1655 + t1627;
t1606 = 0.1e1 / t1612;
t1628 = t1651 * qJ(3,1);
t1613 = t1659 * t1657 + t1628;
t1607 = 0.1e1 / t1613;
t1600 = t1633 * g(1) - t1630 * g(2);
t1673 = g(3) * t1656 + t1600 * t1650;
t1721 = t1673 * t1650;
t1686 = t1649 * t1721;
t1601 = t1634 * g(1) - t1631 * g(2);
t1672 = g(3) * t1658 + t1601 * t1652;
t1723 = t1672 * t1652;
t1690 = t1651 * t1723;
t1599 = t1632 * g(1) - t1629 * g(2);
t1674 = g(3) * t1654 + t1599 * t1648;
t1725 = t1674 * t1648;
t1694 = t1647 * t1725;
t1598 = t1631 * g(1) + t1634 * g(2);
t1774 = -g(3) * t1652 + t1601 * t1658;
t1578 = t1598 * t1651 + t1657 * t1774;
t1662 = 0.1e1 / qJ(3,1);
t1777 = t1578 * t1662;
t1596 = t1629 * g(1) + t1632 * g(2);
t1776 = -g(3) * t1648 + t1599 * t1654;
t1574 = t1596 * t1647 + t1653 * t1776;
t1660 = 0.1e1 / qJ(3,3);
t1778 = t1574 * t1660;
t1597 = t1630 * g(1) + t1633 * g(2);
t1775 = -g(3) * t1650 + t1600 * t1656;
t1570 = t1597 * t1649 + t1655 * t1775;
t1661 = 0.1e1 / qJ(3,2);
t1779 = t1570 * t1661;
t1788 = (t1566 * t1778 - t1629 * t1694) * t1605 + (t1567 * t1779 - t1630 * t1686) * t1606 + (t1568 * t1777 - t1631 * t1690) * t1607;
t1749 = t1629 * qJ(3,3);
t1563 = (t1632 * t1704 - t1749) * t1641 + (t1602 * t1632 + t1629 * t1701) * t1653 + t1749;
t1748 = t1630 * qJ(3,2);
t1564 = (t1633 * t1703 - t1748) * t1642 + (t1603 * t1633 + t1630 * t1700) * t1655 + t1748;
t1747 = t1631 * qJ(3,1);
t1565 = (t1634 * t1702 - t1747) * t1643 + (t1604 * t1634 + t1631 * t1699) * t1657 + t1747;
t1787 = (t1563 * t1778 + t1632 * t1694) * t1605 + (t1564 * t1779 + t1633 * t1686) * t1606 + (t1565 * t1777 + t1634 * t1690) * t1607;
t1766 = t1658 * pkin(4) + t1613 * t1652;
t1717 = t1766 * t1662;
t1681 = t1657 * t1717;
t1768 = t1656 * pkin(4) + t1612 * t1650;
t1718 = t1768 * t1661;
t1682 = t1655 * t1718;
t1770 = t1654 * pkin(4) + t1611 * t1648;
t1719 = t1770 * t1660;
t1683 = t1653 * t1719;
t1780 = (t1574 * t1683 - t1674 * t1707) * t1605 + (t1570 * t1682 - t1673 * t1706) * t1606 + (t1578 * t1681 - t1672 * t1705) * t1607;
t1769 = t1611 * t1654 - t1755;
t1767 = t1612 * t1656 - t1754;
t1765 = t1613 * t1658 - t1753;
t1569 = -t1597 * t1655 + t1649 * t1775;
t1764 = t1569 * t1661;
t1573 = -t1596 * t1653 + t1647 * t1776;
t1763 = t1573 * t1660;
t1577 = -t1598 * t1657 + t1651 * t1774;
t1762 = t1577 * t1662;
t1743 = t1653 * qJ(3,3);
t1742 = t1655 * qJ(3,2);
t1741 = t1657 * qJ(3,1);
t1620 = t1653 * pkin(1) + t1626;
t1560 = -t1596 * t1620 + t1776 * (t1647 * pkin(1) - t1743);
t1740 = t1560 * t1660;
t1621 = t1655 * pkin(1) + t1627;
t1561 = -t1597 * t1621 + t1775 * (t1649 * pkin(1) - t1742);
t1739 = t1561 * t1661;
t1622 = t1657 * pkin(1) + t1628;
t1562 = -t1598 * t1622 + t1774 * (t1651 * pkin(1) - t1741);
t1738 = t1562 * t1662;
t1724 = t1674 * t1654;
t1722 = t1672 * t1658;
t1720 = t1673 * t1656;
t1713 = t1605 * t1648;
t1712 = t1605 * t1654;
t1711 = t1606 * t1650;
t1710 = t1606 * t1656;
t1709 = t1607 * t1652;
t1708 = t1607 * t1658;
t1695 = t1620 * t1725;
t1692 = t1653 * t1725;
t1691 = t1622 * t1723;
t1688 = t1657 * t1723;
t1687 = t1621 * t1721;
t1684 = t1655 * t1721;
t1680 = t1629 * t1713;
t1679 = t1632 * t1713;
t1678 = t1630 * t1711;
t1677 = t1633 * t1711;
t1676 = t1631 * t1709;
t1675 = t1634 * t1709;
t1665 = t1708 * t1774 + t1710 * t1775 + t1712 * t1776;
t1664 = t1676 * t1774 + t1678 * t1775 + t1680 * t1776;
t1663 = t1675 * t1774 + t1677 * t1775 + t1679 * t1776;
t1610 = t1699 - t1741;
t1609 = t1700 - t1742;
t1608 = t1701 - t1743;
t1559 = (-t1577 * t1717 - t1722) * t1657 * t1607 + (-t1569 * t1718 - t1720) * t1655 * t1606 + (-t1573 * t1719 - t1724) * t1653 * t1605;
t1558 = (t1565 * t1762 - t1634 * t1688) * t1607 + (t1564 * t1764 - t1633 * t1684) * t1606 + (t1563 * t1763 - t1632 * t1692) * t1605;
t1557 = (t1568 * t1762 + t1631 * t1688) * t1607 + (t1567 * t1764 + t1630 * t1684) * t1606 + (t1566 * t1763 + t1629 * t1692) * t1605;
t1 = [0, -t1672 * t1675 - t1673 * t1677 - t1674 * t1679, -t1663, 0, 0, 0, 0, 0, t1558, t1787, t1558, t1663, -t1787, -(t1610 * t1631 + t1765 * t1634) * t1762 - (t1609 * t1630 + t1767 * t1633) * t1764 - (t1608 * t1629 + t1769 * t1632) * t1763 + (t1565 * t1738 - t1634 * t1691) * t1607 + (t1564 * t1739 - t1633 * t1687) * t1606 + (t1563 * t1740 - t1632 * t1695) * t1605, -g(1); 0, t1672 * t1676 + t1673 * t1678 + t1674 * t1680, t1664, 0, 0, 0, 0, 0, t1557, t1788, t1557, -t1664, -t1788, -(t1610 * t1634 - t1765 * t1631) * t1762 - (t1609 * t1633 - t1767 * t1630) * t1764 - (t1608 * t1632 - t1769 * t1629) * t1763 + (t1568 * t1738 + t1631 * t1691) * t1607 + (t1567 * t1739 + t1630 * t1687) * t1606 + (t1566 * t1740 + t1629 * t1695) * t1605, -g(2); 0, -t1672 * t1708 - t1673 * t1710 - t1674 * t1712, -t1665, 0, 0, 0, 0, 0, t1559, -t1780, t1559, t1665, t1780, t1766 * t1762 + t1768 * t1764 + t1770 * t1763 + (-t1562 * t1681 - t1622 * t1722) * t1607 + (-t1561 * t1682 - t1621 * t1720) * t1606 + (-t1560 * t1683 - t1620 * t1724) * t1605, -g(3);];
tau_reg  = t1;
