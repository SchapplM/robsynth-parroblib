% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
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
% Datum: 2022-11-04 17:05
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR8V1G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:05:04
% EndTime: 2022-11-04 17:05:05
% DurationCPUTime: 1.25s
% Computational Cost: add. (963->165), mult. (1749->286), div. (159->9), fcn. (1668->26), ass. (0->149)
t1679 = sin(qJ(1,3));
t1685 = cos(qJ(1,3));
t1675 = legFrame(3,2);
t1659 = sin(t1675);
t1662 = cos(t1675);
t1703 = t1662 * g(1) - t1659 * g(2);
t1626 = g(3) * t1679 - t1685 * t1703;
t1672 = pkin(4) + qJ(3,3);
t1668 = 0.1e1 / t1672;
t1738 = t1668 * t1685;
t1724 = t1626 * t1738;
t1652 = cos(pkin(5)) * pkin(2) + pkin(1);
t1678 = sin(qJ(2,3));
t1684 = cos(qJ(2,3));
t1767 = pkin(2) * sin(pkin(5));
t1700 = t1652 * t1684 - t1678 * t1767;
t1632 = 0.1e1 / t1700;
t1753 = t1632 * t1668;
t1778 = t1626 * t1753;
t1681 = sin(qJ(1,2));
t1687 = cos(qJ(1,2));
t1676 = legFrame(2,2);
t1660 = sin(t1676);
t1663 = cos(t1676);
t1702 = t1663 * g(1) - t1660 * g(2);
t1627 = g(3) * t1681 - t1687 * t1702;
t1673 = pkin(4) + qJ(3,2);
t1669 = 0.1e1 / t1673;
t1737 = t1669 * t1687;
t1777 = t1627 * t1737;
t1680 = sin(qJ(2,2));
t1686 = cos(qJ(2,2));
t1699 = t1652 * t1686 - t1680 * t1767;
t1633 = 0.1e1 / t1699;
t1752 = t1633 * t1669;
t1776 = t1627 * t1752;
t1674 = pkin(4) + qJ(3,1);
t1670 = 0.1e1 / t1674;
t1689 = cos(qJ(1,1));
t1736 = t1670 * t1689;
t1683 = sin(qJ(1,1));
t1677 = legFrame(1,2);
t1661 = sin(t1677);
t1664 = cos(t1677);
t1701 = t1664 * g(1) - t1661 * g(2);
t1771 = g(3) * t1683 - t1689 * t1701;
t1775 = t1771 * t1736;
t1682 = sin(qJ(2,1));
t1688 = cos(qJ(2,1));
t1698 = t1652 * t1688 - t1682 * t1767;
t1634 = 0.1e1 / t1698;
t1751 = t1634 * t1670;
t1774 = t1771 * t1751;
t1665 = qJ(2,3) + pkin(5);
t1656 = cos(t1665);
t1770 = pkin(2) * t1656;
t1666 = qJ(2,2) + pkin(5);
t1657 = cos(t1666);
t1769 = pkin(2) * t1657;
t1667 = qJ(2,1) + pkin(5);
t1658 = cos(t1667);
t1768 = pkin(2) * t1658;
t1759 = t1684 * pkin(1);
t1758 = t1686 * pkin(1);
t1757 = t1688 * pkin(1);
t1651 = t1689 * t1757;
t1614 = g(3) * (-t1689 * qJ(3,1) + t1683 * t1757) - t1701 * (t1683 * qJ(3,1) + t1651);
t1756 = t1614 * t1634;
t1649 = t1685 * t1759;
t1615 = g(3) * (-t1685 * qJ(3,3) + t1679 * t1759) - t1703 * (t1679 * qJ(3,3) + t1649);
t1755 = t1615 * t1632;
t1650 = t1687 * t1758;
t1616 = g(3) * (-t1687 * qJ(3,2) + t1681 * t1758) - t1702 * (t1681 * qJ(3,2) + t1650);
t1754 = t1616 * t1633;
t1644 = 0.1e1 / (t1759 + t1770);
t1750 = t1644 * t1659;
t1749 = t1644 * t1662;
t1645 = 0.1e1 / (t1758 + t1769);
t1748 = t1645 * t1660;
t1747 = t1645 * t1663;
t1646 = 0.1e1 / (t1757 + t1768);
t1746 = t1646 * t1661;
t1745 = t1646 * t1664;
t1744 = t1652 * t1662;
t1743 = t1652 * t1663;
t1742 = t1652 * t1664;
t1741 = t1659 * t1652;
t1740 = t1660 * t1652;
t1739 = t1661 * t1652;
t1735 = t1662 * t1767;
t1734 = t1663 * t1767;
t1733 = t1664 * t1767;
t1732 = t1659 * t1767;
t1731 = t1660 * t1767;
t1730 = t1661 * t1767;
t1697 = g(3) * t1685 + t1679 * t1703;
t1720 = t1697 * t1753;
t1696 = g(3) * t1687 + t1681 * t1702;
t1719 = t1696 * t1752;
t1695 = g(3) * t1689 + t1683 * t1701;
t1718 = t1695 * t1751;
t1715 = t1658 * t1774;
t1714 = t1688 * t1774;
t1713 = t1657 * t1776;
t1712 = t1680 * t1776;
t1711 = t1686 * t1776;
t1653 = sin(t1665);
t1710 = t1653 * t1778;
t1709 = t1656 * t1778;
t1708 = t1678 * t1778;
t1707 = t1684 * t1778;
t1654 = sin(t1666);
t1706 = t1654 * t1776;
t1655 = sin(t1667);
t1705 = t1655 * t1774;
t1704 = t1682 * t1774;
t1641 = t1659 * g(1) + t1662 * g(2);
t1608 = -t1641 * t1684 + t1678 * t1697;
t1642 = t1660 * g(1) + t1663 * g(2);
t1610 = -t1642 * t1686 + t1680 * t1696;
t1643 = t1661 * g(1) + t1664 * g(2);
t1612 = -t1643 * t1688 + t1682 * t1695;
t1694 = t1608 * t1750 + t1610 * t1748 + t1612 * t1746;
t1693 = t1608 * t1749 + t1610 * t1747 + t1612 * t1745;
t1692 = t1695 * t1736 + t1696 * t1737 + t1697 * t1738;
t1596 = (-t1679 * t1741 + t1735) * t1684 + (t1679 * t1732 + t1744) * t1678;
t1597 = (-t1681 * t1740 + t1734) * t1686 + (t1681 * t1731 + t1743) * t1680;
t1598 = (-t1683 * t1739 + t1733) * t1688 + (t1683 * t1730 + t1742) * t1682;
t1691 = t1596 * t1720 + t1597 * t1719 + t1598 * t1718;
t1599 = (t1679 * t1744 + t1732) * t1684 + (-t1679 * t1735 + t1741) * t1678;
t1600 = (t1681 * t1743 + t1731) * t1686 + (-t1681 * t1734 + t1740) * t1680;
t1601 = (t1683 * t1742 + t1730) * t1688 + (-t1683 * t1733 + t1739) * t1682;
t1690 = t1599 * t1720 + t1600 * t1719 + t1601 * t1718;
t1637 = t1682 * t1652 + t1688 * t1767;
t1636 = t1680 * t1652 + t1686 * t1767;
t1635 = t1678 * t1652 + t1684 * t1767;
t1619 = -t1689 * t1674 + t1683 * t1698;
t1618 = -t1687 * t1673 + t1681 * t1699;
t1617 = -t1685 * t1672 + t1679 * t1700;
t1613 = t1643 * t1682 + t1688 * t1695;
t1611 = t1642 * t1680 + t1686 * t1696;
t1609 = t1641 * t1678 + t1684 * t1697;
t1607 = t1643 * t1655 + t1658 * t1695;
t1606 = -t1643 * t1658 + t1655 * t1695;
t1605 = t1642 * t1654 + t1657 * t1696;
t1604 = -t1642 * t1657 + t1654 * t1696;
t1603 = t1641 * t1653 + t1656 * t1697;
t1602 = -t1641 * t1656 + t1653 * t1697;
t1 = [0, t1599 * t1778 + t1600 * t1776 + t1601 * t1774, t1690, 0, 0, 0, 0, 0, t1599 * t1707 + t1600 * t1711 + t1601 * t1714 + t1694, -t1599 * t1708 - t1600 * t1712 - t1601 * t1704 + t1609 * t1750 + t1611 * t1748 + t1613 * t1746, t1599 * t1709 + t1600 * t1713 + t1601 * t1715 + t1602 * t1750 + t1604 * t1748 + t1606 * t1746, -t1599 * t1710 - t1600 * t1706 - t1601 * t1705 + t1603 * t1750 + t1605 * t1748 + t1607 * t1746, -t1690, (t1601 * t1756 - (t1619 * t1664 + t1661 * t1637) * t1771) * t1670 + (t1600 * t1754 - (t1618 * t1663 + t1660 * t1636) * t1627) * t1669 + (t1599 * t1755 - (t1617 * t1662 + t1659 * t1635) * t1626) * t1668 + t1694 * pkin(1), -g(1); 0, t1596 * t1778 + t1597 * t1776 + t1598 * t1774, t1691, 0, 0, 0, 0, 0, t1596 * t1707 + t1597 * t1711 + t1598 * t1714 + t1693, -t1596 * t1708 - t1597 * t1712 - t1598 * t1704 + t1609 * t1749 + t1611 * t1747 + t1613 * t1745, t1596 * t1709 + t1597 * t1713 + t1598 * t1715 + t1602 * t1749 + t1604 * t1747 + t1606 * t1745, -t1596 * t1710 - t1597 * t1706 - t1598 * t1705 + t1603 * t1749 + t1605 * t1747 + t1607 * t1745, -t1691, (t1598 * t1756 - (-t1619 * t1661 + t1664 * t1637) * t1771) * t1670 + (t1597 * t1754 - (-t1618 * t1660 + t1663 * t1636) * t1627) * t1669 + (t1596 * t1755 - (-t1617 * t1659 + t1662 * t1635) * t1626) * t1668 + t1693 * pkin(1), -g(2); 0, t1777 + t1775 + t1724, t1692, 0, 0, 0, 0, 0, t1684 * t1724 + t1686 * t1777 + t1688 * t1775, -t1678 * t1724 - t1680 * t1777 - t1682 * t1775, t1656 * t1724 + t1657 * t1777 + t1658 * t1775, -t1653 * t1724 - t1654 * t1777 - t1655 * t1775, -t1692, (t1689 * t1614 - (t1683 * t1674 + t1689 * t1768 + t1651) * t1771) * t1670 + (t1687 * t1616 - (t1681 * t1673 + t1687 * t1769 + t1650) * t1627) * t1669 + (t1685 * t1615 - (t1679 * t1672 + t1685 * t1770 + t1649) * t1626) * t1668, -g(3);];
tau_reg  = t1;
