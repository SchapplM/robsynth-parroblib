% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [13x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR1V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:34
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR1V1G2A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(13,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_mdp: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_mdp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:34:15
% EndTime: 2020-08-06 19:34:18
% DurationCPUTime: 2.81s
% Computational Cost: add. (7050->243), mult. (11332->537), div. (1548->21), fcn. (8961->18), ass. (0->229)
t1830 = 2 * pkin(1);
t1690 = pkin(1) ^ 2;
t1829 = -2 * t1690;
t1828 = 2 * MDP(5);
t1668 = pkin(3) + qJ(3,3);
t1637 = 0.1e1 / t1668;
t1669 = pkin(3) + qJ(3,2);
t1640 = 0.1e1 / t1669;
t1670 = pkin(3) + qJ(3,1);
t1643 = 0.1e1 / t1670;
t1680 = cos(qJ(2,3));
t1650 = 0.1e1 / t1680;
t1682 = cos(qJ(2,2));
t1655 = 0.1e1 / t1682;
t1684 = cos(qJ(2,1));
t1660 = 0.1e1 / t1684;
t1689 = pkin(1) + pkin(2);
t1665 = 1 / t1689;
t1638 = 0.1e1 / t1668 ^ 2;
t1641 = 0.1e1 / t1669 ^ 2;
t1644 = 0.1e1 / t1670 ^ 2;
t1666 = 1 / t1689 ^ 2;
t1827 = MDP(6) * t1666;
t1664 = t1689 ^ 2;
t1826 = MDP(8) * t1665 / t1664;
t1679 = sin(qJ(1,1));
t1685 = cos(qJ(1,1));
t1775 = t1684 * t1689;
t1627 = t1670 * t1679 + t1685 * t1775;
t1673 = legFrame(1,2);
t1633 = sin(t1673);
t1636 = cos(t1673);
t1678 = sin(qJ(2,1));
t1686 = xDP(3);
t1687 = xDP(2);
t1688 = xDP(1);
t1778 = t1679 * t1684;
t1787 = t1670 * t1685;
t1600 = (-t1688 * t1787 + (t1678 * t1687 + t1688 * t1778) * t1689) * t1636 + (t1687 * t1787 - (-t1678 * t1688 + t1687 * t1778) * t1689) * t1633 + t1627 * t1686;
t1624 = t1633 * t1688 + t1636 * t1687;
t1780 = t1678 * t1624;
t1606 = (t1685 * t1686 + (-t1633 * t1687 + t1636 * t1688) * t1679) * t1684 + t1780;
t1621 = t1624 ^ 2;
t1645 = t1643 * t1644;
t1659 = t1684 ^ 2;
t1662 = t1660 / t1659;
t1796 = t1644 * t1660;
t1576 = -t1606 * t1645 * t1660 * t1600 + (-(-t1606 * t1689 + t1600) * t1606 * t1796 + t1665 * t1621 * t1662) * t1643;
t1825 = qJ(3,1) * t1576;
t1677 = sin(qJ(1,2));
t1683 = cos(qJ(1,2));
t1776 = t1682 * t1689;
t1626 = t1669 * t1677 + t1683 * t1776;
t1672 = legFrame(2,2);
t1632 = sin(t1672);
t1635 = cos(t1672);
t1676 = sin(qJ(2,2));
t1781 = t1677 * t1682;
t1788 = t1669 * t1683;
t1599 = (-t1688 * t1788 + (t1676 * t1687 + t1688 * t1781) * t1689) * t1635 + (t1687 * t1788 - (-t1676 * t1688 + t1687 * t1781) * t1689) * t1632 + t1626 * t1686;
t1623 = t1632 * t1688 + t1635 * t1687;
t1783 = t1676 * t1623;
t1605 = (t1683 * t1686 + (-t1632 * t1687 + t1635 * t1688) * t1677) * t1682 + t1783;
t1620 = t1623 ^ 2;
t1642 = t1640 * t1641;
t1654 = t1682 ^ 2;
t1657 = t1655 / t1654;
t1797 = t1641 * t1655;
t1575 = -t1605 * t1642 * t1655 * t1599 + (-(-t1605 * t1689 + t1599) * t1605 * t1797 + t1665 * t1620 * t1657) * t1640;
t1824 = qJ(3,2) * t1575;
t1675 = sin(qJ(1,3));
t1681 = cos(qJ(1,3));
t1777 = t1680 * t1689;
t1625 = t1668 * t1675 + t1681 * t1777;
t1671 = legFrame(3,2);
t1631 = sin(t1671);
t1634 = cos(t1671);
t1674 = sin(qJ(2,3));
t1784 = t1675 * t1680;
t1789 = t1668 * t1681;
t1598 = (-t1688 * t1789 + (t1674 * t1687 + t1688 * t1784) * t1689) * t1634 + (t1687 * t1789 - (-t1674 * t1688 + t1687 * t1784) * t1689) * t1631 + t1625 * t1686;
t1622 = t1631 * t1688 + t1634 * t1687;
t1786 = t1674 * t1622;
t1604 = (t1681 * t1686 + (-t1631 * t1687 + t1634 * t1688) * t1675) * t1680 + t1786;
t1619 = t1622 ^ 2;
t1639 = t1637 * t1638;
t1649 = t1680 ^ 2;
t1652 = t1650 / t1649;
t1798 = t1638 * t1650;
t1574 = -t1604 * t1639 * t1650 * t1598 + (-(-t1604 * t1689 + t1598) * t1604 * t1798 + t1665 * t1619 * t1652) * t1637;
t1823 = qJ(3,3) * t1574;
t1822 = t1598 * t1689;
t1821 = t1599 * t1689;
t1820 = t1600 * t1689;
t1601 = t1604 ^ 2;
t1819 = t1601 * t1638;
t1602 = t1605 ^ 2;
t1818 = t1602 * t1641;
t1603 = t1606 ^ 2;
t1817 = t1603 * t1644;
t1816 = t1604 * t1638;
t1815 = t1604 * t1674;
t1814 = t1605 * t1641;
t1813 = t1605 * t1676;
t1812 = t1606 * t1644;
t1811 = t1606 * t1678;
t1810 = t1619 * t1637;
t1809 = t1619 * t1666;
t1808 = t1620 * t1640;
t1807 = t1620 * t1666;
t1806 = t1621 * t1643;
t1805 = t1621 * t1666;
t1804 = t1631 * t1650;
t1803 = t1632 * t1655;
t1802 = t1633 * t1660;
t1801 = t1634 * t1650;
t1800 = t1635 * t1655;
t1799 = t1636 * t1660;
t1795 = t1650 * t1674;
t1651 = 0.1e1 / t1680 ^ 2;
t1794 = t1651 * t1674;
t1793 = t1655 * t1676;
t1656 = 0.1e1 / t1682 ^ 2;
t1792 = t1656 * t1676;
t1791 = t1660 * t1678;
t1661 = 0.1e1 / t1684 ^ 2;
t1790 = t1661 * t1678;
t1785 = t1674 * t1689;
t1782 = t1676 * t1689;
t1779 = t1678 * t1689;
t1774 = t1674 * t1828;
t1773 = t1676 * t1828;
t1772 = t1678 * t1828;
t1771 = t1574 * t1795;
t1770 = t1575 * t1793;
t1769 = t1576 * t1791;
t1768 = t1601 * t1639 * t1651;
t1767 = t1602 * t1642 * t1656;
t1766 = t1603 * t1645 * t1661;
t1765 = t1622 * t1816;
t1764 = t1637 * t1815;
t1763 = t1623 * t1814;
t1762 = t1640 * t1813;
t1761 = t1624 * t1812;
t1760 = t1643 * t1811;
t1759 = t1652 * t1809;
t1653 = 0.1e1 / t1649 ^ 2;
t1758 = t1619 * t1653 * t1674;
t1757 = t1657 * t1807;
t1658 = 0.1e1 / t1654 ^ 2;
t1756 = t1620 * t1658 * t1676;
t1755 = t1662 * t1805;
t1663 = 0.1e1 / t1659 ^ 2;
t1754 = t1621 * t1663 * t1678;
t1646 = t1674 ^ 2;
t1753 = MDP(4) * t1646 + MDP(1);
t1647 = t1676 ^ 2;
t1752 = MDP(4) * t1647 + MDP(1);
t1648 = t1678 ^ 2;
t1751 = MDP(4) * t1648 + MDP(1);
t1711 = (t1622 - t1815) * t1650 ^ 2 * t1637 * t1622 + (-((-t1604 + t1786) * t1668 * t1650 + (-t1604 * t1664 + t1822) * t1680 * t1637) * t1798 - t1639 * t1822) * t1604;
t1732 = t1622 * t1665 * t1764;
t1750 = (-pkin(1) * t1574 * t1680 + (-qJ(3,3) * t1819 + t1732 * t1830) * t1651 + t1711) * MDP(12);
t1710 = (t1623 - t1813) * t1655 ^ 2 * t1640 * t1623 + (-((-t1605 + t1783) * t1669 * t1655 + (-t1605 * t1664 + t1821) * t1682 * t1640) * t1797 - t1642 * t1821) * t1605;
t1731 = t1623 * t1665 * t1762;
t1749 = (-pkin(1) * t1575 * t1682 + (-qJ(3,2) * t1818 + t1731 * t1830) * t1656 + t1710) * MDP(12);
t1709 = (t1624 - t1811) * t1660 ^ 2 * t1643 * t1624 + (-((-t1606 + t1780) * t1670 * t1660 + (-t1606 * t1664 + t1820) * t1684 * t1643) * t1796 - t1645 * t1820) * t1606;
t1730 = t1624 * t1665 * t1760;
t1748 = (-pkin(1) * t1576 * t1684 + (-qJ(3,1) * t1817 + t1730 * t1830) * t1661 + t1709) * MDP(12);
t1747 = 0.2e1 * t1765;
t1746 = 0.2e1 * t1763;
t1745 = 0.2e1 * t1761;
t1744 = t1794 * t1819;
t1743 = t1792 * t1818;
t1742 = t1790 * t1817;
t1741 = t1681 * t1765;
t1740 = t1683 * t1763;
t1739 = t1685 * t1761;
t1738 = t1646 * t1759;
t1737 = t1647 * t1757;
t1736 = t1648 * t1755;
t1735 = (pkin(1) * t1604 - 0.2e1 * t1598) * t1637 * t1651 * t1764;
t1734 = (pkin(1) * t1605 - 0.2e1 * t1599) * t1640 * t1656 * t1762;
t1733 = (pkin(1) * t1606 - 0.2e1 * t1600) * t1643 * t1661 * t1760;
t1628 = 0.2e1 * t1649 - 0.1e1;
t1729 = t1628 * t1652 * t1747;
t1728 = t1747 * t1794;
t1629 = 0.2e1 * t1654 - 0.1e1;
t1727 = t1629 * t1657 * t1746;
t1726 = t1746 * t1792;
t1630 = 0.2e1 * t1659 - 0.1e1;
t1725 = t1630 * t1662 * t1745;
t1724 = t1745 * t1790;
t1723 = (t1646 * t1653 + t1651) * t1810;
t1722 = (t1647 * t1658 + t1656) * t1808;
t1721 = (t1648 * t1663 + t1661) * t1806;
t1720 = t1675 * t1777 - t1789;
t1719 = t1677 * t1776 - t1788;
t1718 = t1679 * t1775 - t1787;
t1717 = -pkin(1) * t1809 + 0.2e1 * t1598 * t1816;
t1716 = -pkin(1) * t1807 + 0.2e1 * t1599 * t1814;
t1715 = -pkin(1) * t1805 + 0.2e1 * t1600 * t1812;
t1714 = ((qJ(3,3) ^ 2 + t1649 * t1690) * t1574 + (-qJ(3,3) * t1738 - t1680 * t1711) * pkin(1) + (t1717 * qJ(3,3) + t1732 * t1829) * t1650) * MDP(12) + (-pkin(1) * t1738 + t1717 * t1650 + 0.2e1 * t1823) * MDP(11);
t1713 = ((qJ(3,2) ^ 2 + t1654 * t1690) * t1575 + (-qJ(3,2) * t1737 - t1682 * t1710) * pkin(1) + (t1716 * qJ(3,2) + t1731 * t1829) * t1655) * MDP(12) + (-pkin(1) * t1737 + t1716 * t1655 + 0.2e1 * t1824) * MDP(11);
t1712 = ((qJ(3,1) ^ 2 + t1659 * t1690) * t1576 + (-qJ(3,1) * t1736 - t1684 * t1709) * pkin(1) + (t1715 * qJ(3,1) + t1730 * t1829) * t1660) * MDP(12) + (-pkin(1) * t1736 + t1715 * t1660 + 0.2e1 * t1825) * MDP(11);
t1708 = t1631 * t1771 + t1632 * t1770 + t1633 * t1769;
t1707 = t1634 * t1771 + t1635 * t1770 + t1636 * t1769;
t1706 = t1574 * t1774 + (t1753 * t1574 + t1714) * t1650;
t1705 = t1575 * t1773 + (t1752 * t1575 + t1713) * t1655;
t1704 = t1576 * t1772 + (t1751 * t1576 + t1712) * t1660;
t1618 = t1633 * t1678 + t1636 * t1778;
t1617 = t1632 * t1676 + t1635 * t1781;
t1616 = t1631 * t1674 + t1634 * t1784;
t1615 = -t1633 * t1778 + t1636 * t1678;
t1614 = -t1632 * t1781 + t1635 * t1676;
t1613 = -t1631 * t1784 + t1634 * t1674;
t1612 = t1633 * t1779 + t1718 * t1636;
t1611 = t1632 * t1782 + t1719 * t1635;
t1610 = t1631 * t1785 + t1720 * t1634;
t1609 = -t1719 * t1632 + t1635 * t1782;
t1608 = -t1718 * t1633 + t1636 * t1779;
t1607 = -t1720 * t1631 + t1634 * t1785;
t1594 = (t1661 - 0.2e1) * t1817;
t1593 = (t1656 - 0.2e1) * t1818;
t1592 = (t1651 - 0.2e1) * t1819;
t1570 = (pkin(1) * t1755 - t1825) * t1678;
t1569 = (pkin(1) * t1757 - t1824) * t1676;
t1568 = (pkin(1) * t1759 - t1823) * t1674;
t1 = [(-t1610 * t1768 - t1611 * t1767 - t1612 * t1766) * MDP(11) + (t1631 * t1758 + t1632 * t1756 + t1633 * t1754) * t1826 + (t1612 * t1748 + t1704 * t1618) * t1643 + (t1611 * t1749 + t1705 * t1617) * t1640 + (t1610 * t1750 + t1706 * t1616) * t1637 + (t1616 * t1723 + t1617 * t1722 + t1618 * t1721) * t1827 + ((t1616 * t1728 + t1617 * t1726 + t1618 * t1724 - t1631 * t1744 - t1632 * t1743 - t1633 * t1742) * MDP(4) + (t1592 * t1804 + t1593 * t1803 + t1594 * t1802 + t1616 * t1729 + t1617 * t1727 + t1618 * t1725) * MDP(5) + t1708 * MDP(6) + (t1574 * t1631 + t1575 * t1632 + t1576 * t1633) * MDP(7) + (-t1708 * MDP(11) + (t1568 * t1804 + t1569 * t1803 + t1570 * t1802 + t1631 * t1735 + t1632 * t1734 + t1633 * t1733) * MDP(12)) * pkin(1)) * t1665; (-t1607 * t1768 - t1608 * t1766 - t1609 * t1767) * MDP(11) + (t1634 * t1758 + t1635 * t1756 + t1636 * t1754) * t1826 + (t1608 * t1748 + t1704 * t1615) * t1643 + (t1609 * t1749 + t1705 * t1614) * t1640 + (t1607 * t1750 + t1706 * t1613) * t1637 + (t1613 * t1723 + t1614 * t1722 + t1615 * t1721) * t1827 + ((t1613 * t1728 + t1614 * t1726 + t1615 * t1724 - t1634 * t1744 - t1635 * t1743 - t1636 * t1742) * MDP(4) + (t1592 * t1801 + t1593 * t1800 + t1594 * t1799 + t1613 * t1729 + t1614 * t1727 + t1615 * t1725) * MDP(5) + t1707 * MDP(6) + (t1574 * t1634 + t1575 * t1635 + t1576 * t1636) * MDP(7) + (-t1707 * MDP(11) + (t1568 * t1801 + t1569 * t1800 + t1570 * t1799 + t1634 * t1735 + t1635 * t1734 + t1636 * t1733) * MDP(12)) * pkin(1)) * t1665; (-t1625 * t1768 - t1626 * t1767 - t1627 * t1766) * MDP(11) + (0.2e1 * (t1739 * t1791 + t1740 * t1793 + t1741 * t1795) * MDP(4) + (t1628 * t1651 * t1741 + t1629 * t1656 * t1740 + t1630 * t1661 * t1739) * t1828) * t1665 + (t1627 * t1748 + ((t1684 * t1772 + t1751) * t1576 + t1712) * t1685) * t1643 + (t1626 * t1749 + ((t1682 * t1773 + t1752) * t1575 + t1713) * t1683) * t1640 + (t1625 * t1750 + ((t1680 * t1774 + t1753) * t1574 + t1714) * t1681) * t1637 + ((t1648 * t1662 + t1660) * t1685 * t1806 + (t1647 * t1657 + t1655) * t1683 * t1808 + (t1646 * t1652 + t1650) * t1681 * t1810) * t1827;];
taucX  = t1;
