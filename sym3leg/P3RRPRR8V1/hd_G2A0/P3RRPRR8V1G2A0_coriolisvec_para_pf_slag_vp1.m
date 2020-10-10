% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR8V1G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:59:05
% EndTime: 2020-08-06 19:59:08
% DurationCPUTime: 2.79s
% Computational Cost: add. (13494->288), mult. (22866->474), div. (2937->12), fcn. (18921->44), ass. (0->233)
t730 = pkin(4) + qJ(3,3);
t737 = sin(qJ(1,3));
t743 = cos(qJ(1,3));
t726 = cos(pkin(5));
t828 = pkin(2) * t726;
t676 = pkin(1) + t828;
t742 = cos(qJ(2,3));
t725 = sin(pkin(5));
t736 = sin(qJ(2,3));
t815 = t725 * t736;
t768 = pkin(2) * t815 - t676 * t742;
t633 = -t730 * t743 - t768 * t737;
t829 = pkin(2) * t725;
t645 = t676 * t736 + t742 * t829;
t733 = legFrame(3,2);
t695 = sin(t733);
t698 = cos(t733);
t611 = t633 * t698 + t645 * t695;
t612 = -t633 * t695 + t645 * t698;
t713 = qJ(2,3) + pkin(5);
t692 = cos(t713);
t705 = t742 * pkin(1);
t851 = pkin(2) * t692 + t705;
t636 = t730 * t737 + t743 * t851;
t716 = 0.1e1 / t730;
t753 = xDP(3);
t754 = xDP(2);
t755 = xDP(1);
t599 = (t611 * t755 + t612 * t754 + t636 * t753) * t716;
t864 = 0.2e1 * t599;
t731 = pkin(4) + qJ(3,2);
t739 = sin(qJ(1,2));
t745 = cos(qJ(1,2));
t744 = cos(qJ(2,2));
t738 = sin(qJ(2,2));
t814 = t725 * t738;
t767 = pkin(2) * t814 - t676 * t744;
t634 = -t731 * t745 - t767 * t739;
t646 = t676 * t738 + t744 * t829;
t734 = legFrame(2,2);
t696 = sin(t734);
t699 = cos(t734);
t613 = t634 * t699 + t646 * t696;
t614 = -t634 * t696 + t646 * t699;
t714 = qJ(2,2) + pkin(5);
t693 = cos(t714);
t706 = t744 * pkin(1);
t850 = pkin(2) * t693 + t706;
t637 = t731 * t739 + t745 * t850;
t717 = 0.1e1 / t731;
t600 = (t613 * t755 + t614 * t754 + t637 * t753) * t717;
t863 = 0.2e1 * t600;
t732 = pkin(4) + qJ(3,1);
t741 = sin(qJ(1,1));
t747 = cos(qJ(1,1));
t746 = cos(qJ(2,1));
t740 = sin(qJ(2,1));
t813 = t725 * t740;
t766 = pkin(2) * t813 - t676 * t746;
t635 = -t732 * t747 - t766 * t741;
t647 = t676 * t740 + t746 * t829;
t735 = legFrame(1,2);
t697 = sin(t735);
t700 = cos(t735);
t615 = t635 * t700 + t647 * t697;
t616 = -t635 * t697 + t647 * t700;
t715 = qJ(2,1) + pkin(5);
t694 = cos(t715);
t707 = t746 * pkin(1);
t849 = pkin(2) * t694 + t707;
t638 = t732 * t741 + t747 * t849;
t718 = 0.1e1 / t732;
t601 = (t615 * t755 + t616 * t754 + t638 * t753) * t718;
t862 = 0.2e1 * t601;
t803 = 0.2e1 * pkin(1);
t648 = 0.1e1 / t851;
t626 = (t695 * t755 + t698 * t754) * t648;
t623 = t626 ^ 2;
t763 = pkin(2) ^ 2;
t764 = pkin(1) ^ 2;
t804 = -t763 - t764;
t658 = t803 * t828 - t804;
t860 = t623 * t658;
t649 = 0.1e1 / t850;
t627 = (t696 * t755 + t699 * t754) * t649;
t624 = t627 ^ 2;
t859 = t624 * t658;
t650 = 0.1e1 / t849;
t628 = (t697 * t755 + t700 * t754) * t650;
t625 = t628 ^ 2;
t858 = t625 * t658;
t686 = sin(t713);
t702 = t736 * pkin(1);
t857 = rSges(3,1) * t686 + rSges(3,2) * t692 + t702;
t856 = pkin(2) * t686 + t702;
t687 = sin(t714);
t703 = t738 * pkin(1);
t855 = rSges(3,1) * t687 + rSges(3,2) * t693 + t703;
t854 = pkin(2) * t687 + t703;
t688 = sin(t715);
t704 = t740 * pkin(1);
t853 = rSges(3,1) * t688 + rSges(3,2) * t694 + t704;
t852 = pkin(2) * t688 + t704;
t848 = -2 * m(3);
t847 = 2 * m(3);
t846 = pkin(1) * m(3);
t845 = m(2) * rSges(2,3);
t844 = m(3) * rSges(3,2);
t760 = rSges(2,2) ^ 2;
t762 = rSges(2,1) ^ 2;
t651 = m(3) * t764 + (-t760 + t762) * m(2) + Icges(2,2) - Icges(2,1);
t843 = -t651 / 0.2e1;
t759 = rSges(3,2) ^ 2;
t761 = rSges(3,1) ^ 2;
t660 = m(3) * (-t759 + t761) - Icges(3,1) + Icges(3,2);
t842 = -t660 / 0.2e1;
t841 = t660 / 0.2e1;
t680 = -rSges(3,1) * t844 + Icges(3,4);
t840 = -t680 / 0.2e1;
t682 = m(2) * rSges(2,1) * rSges(2,2) - Icges(2,4);
t839 = t682 / 0.2e1;
t727 = (qJ(3,3) + rSges(3,3));
t838 = m(3) * t727;
t728 = (qJ(3,2) + rSges(3,3));
t837 = m(3) * t728;
t729 = (qJ(3,1) + rSges(3,3));
t836 = m(3) * t729;
t756 = 0.2e1 * qJ(2,3);
t719 = sin(t756);
t827 = t651 * t719;
t757 = 0.2e1 * qJ(2,2);
t720 = sin(t757);
t826 = t651 * t720;
t758 = 0.2e1 * qJ(2,1);
t721 = sin(t758);
t825 = t651 * t721;
t824 = t676 * t737;
t823 = t676 * t739;
t822 = t676 * t741;
t677 = 0.2e1 * t713;
t666 = cos(t677);
t821 = t680 * t666;
t678 = 0.2e1 * t714;
t667 = cos(t678);
t820 = t680 * t667;
t679 = 0.2e1 * t715;
t668 = cos(t679);
t819 = t680 * t668;
t722 = cos(t756);
t818 = t682 * t722;
t723 = cos(t757);
t817 = t682 * t723;
t724 = cos(t758);
t816 = t682 * t724;
t812 = -rSges(2,1) * t845 + Icges(2,5);
t681 = rSges(2,2) * t845 - Icges(2,6);
t801 = t737 * t829;
t617 = (-t695 * t824 + t698 * t829) * t742 + t736 * (t676 * t698 + t695 * t801);
t620 = (t695 * t829 + t698 * t824) * t742 + (t676 * t695 - t698 * t801) * t736;
t642 = 0.1e1 / t768;
t596 = (t743 * t753 - (t617 * t754 + t620 * t755) * t642) * t716;
t797 = pkin(2) * t803;
t710 = pkin(5) + t756;
t689 = cos(t710);
t808 = -t689 - t726;
t578 = (-t860 + ((-t763 * t666 - t722 * t764 + t808 * t797 + t804) * t596 / 0.2e1 + t851 * t864 + (-t596 * t730 + 0.2e1 * t626 * t856) * t730) * t596) * t716;
t590 = (-0.1e1 / (t705 + (t726 * t742 - t815) * pkin(2)) * t860 + (t768 * t596 + t864) * t596) * t716;
t669 = rSges(3,2) * t838 - Icges(3,6);
t672 = rSges(3,1) * t838 - Icges(3,5);
t774 = pkin(1) * t838 - t812;
t602 = -(t669 * t725 - t672 * t726 - t774) * t736 + t742 * (t669 * t726 + t672 * t725 + t681);
t639 = rSges(3,1) * t692 - rSges(3,2) * t686 + t705;
t663 = sin(t677);
t683 = sin(t710);
t675 = t725 * pkin(1) * t844;
t805 = t760 + t762;
t765 = -(rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - ((2 * rSges(2,3) ^ 2) + t805) * m(2) / 0.2e1 + t675 - Icges(3,2) / 0.2e1 - Icges(2,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,1) / 0.2e1 - Icges(1,3);
t780 = rSges(3,1) * t683 + rSges(3,2) * t689;
t781 = -t759 / 0.2e1 - t761 / 0.2e1 - t764 / 0.2e1;
t796 = t623 * t648 * t856;
t798 = m(3) * t803;
t811 = (t666 * t842 + t722 * t843 + t765) * t590 - t602 * t796 + 0.2e1 * (t663 * t840 + t719 * t839) * t590 + (t639 * t578 + (-(t727 ^ 2) + t781 + (t808 * rSges(3,1) + rSges(3,2) * t683) * pkin(1)) * t590) * m(3) + (t669 * t686 - t672 * t692 + t681 * t736 - t774 * t742) * t623 + (t838 * t864 + (-t660 * t663 - t780 * t798 - 0.2e1 * t818 + 0.2e1 * t821 - t827) * t626) * t596;
t800 = t739 * t829;
t618 = (-t696 * t823 + t699 * t829) * t744 + t738 * (t676 * t699 + t696 * t800);
t621 = (t696 * t829 + t699 * t823) * t744 + (t676 * t696 - t699 * t800) * t738;
t643 = 0.1e1 / t767;
t597 = (t745 * t753 - (t618 * t754 + t621 * t755) * t643) * t717;
t711 = pkin(5) + t757;
t690 = cos(t711);
t807 = -t690 - t726;
t579 = (-t859 + ((-t763 * t667 - t723 * t764 + t807 * t797 + t804) * t597 / 0.2e1 + t850 * t863 + (-t597 * t731 + 0.2e1 * t627 * t854) * t731) * t597) * t717;
t591 = (-0.1e1 / (t706 + (t726 * t744 - t814) * pkin(2)) * t859 + (t767 * t597 + t863) * t597) * t717;
t670 = rSges(3,2) * t837 - Icges(3,6);
t673 = rSges(3,1) * t837 - Icges(3,5);
t773 = pkin(1) * t837 - t812;
t603 = -(t670 * t725 - t673 * t726 - t773) * t738 + t744 * (t670 * t726 + t673 * t725 + t681);
t640 = rSges(3,1) * t693 - rSges(3,2) * t687 + t706;
t664 = sin(t678);
t684 = sin(t711);
t779 = rSges(3,1) * t684 + rSges(3,2) * t690;
t795 = t624 * t649 * t854;
t810 = (t667 * t842 + t723 * t843 + t765) * t591 - t603 * t795 + 0.2e1 * (t664 * t840 + t720 * t839) * t591 + (t640 * t579 + (-(t728 ^ 2) + t781 + (t807 * rSges(3,1) + rSges(3,2) * t684) * pkin(1)) * t591) * m(3) + (t670 * t687 - t673 * t693 + t681 * t738 - t773 * t744) * t624 + (t837 * t863 + (-t660 * t664 - t779 * t798 - 0.2e1 * t817 + 0.2e1 * t820 - t826) * t627) * t597;
t799 = t741 * t829;
t619 = (-t697 * t822 + t700 * t829) * t746 + t740 * (t676 * t700 + t697 * t799);
t622 = (t697 * t829 + t700 * t822) * t746 + (t676 * t697 - t700 * t799) * t740;
t644 = 0.1e1 / t766;
t598 = (t747 * t753 - (t619 * t754 + t622 * t755) * t644) * t718;
t712 = pkin(5) + t758;
t691 = cos(t712);
t806 = -t691 - t726;
t580 = (-t858 + ((-t763 * t668 - t724 * t764 + t806 * t797 + t804) * t598 / 0.2e1 + t849 * t862 + (-t598 * t732 + 0.2e1 * t628 * t852) * t732) * t598) * t718;
t592 = (-0.1e1 / (t707 + (t726 * t746 - t813) * pkin(2)) * t858 + (t766 * t598 + t862) * t598) * t718;
t671 = rSges(3,2) * t836 - Icges(3,6);
t674 = rSges(3,1) * t836 - Icges(3,5);
t772 = pkin(1) * t836 - t812;
t604 = -(t671 * t725 - t674 * t726 - t772) * t740 + t746 * (t671 * t726 + t674 * t725 + t681);
t641 = rSges(3,1) * t694 - rSges(3,2) * t688 + t707;
t665 = sin(t679);
t685 = sin(t712);
t778 = rSges(3,1) * t685 + rSges(3,2) * t691;
t794 = t625 * t650 * t852;
t809 = (t668 * t842 + t724 * t843 + t765) * t592 - t604 * t794 + 0.2e1 * (t665 * t840 + t721 * t839) * t592 + (t641 * t580 + (-(t729 ^ 2) + t781 + (t806 * rSges(3,1) + rSges(3,2) * t685) * pkin(1)) * t592) * m(3) + (t671 * t688 - t674 * t694 + t681 * t740 - t772 * t746) * t625 + (t836 * t862 + (-t660 * t665 - t778 * t798 - 0.2e1 * t816 + 0.2e1 * t819 - t825) * t628) * t598;
t793 = t811 * t642;
t792 = t810 * t643;
t791 = t809 * t644;
t790 = (-t596 * t727 / 0.2e1 + t857 * t626) * t596 * t847 + (t590 * t639 - t578) * m(3);
t789 = (-t597 * t728 / 0.2e1 + t855 * t627) * t597 * t847 + (t591 * t640 - t579) * m(3);
t788 = (-t598 * t729 / 0.2e1 + t853 * t628) * t598 * t847 + (t592 * t641 - t580) * m(3);
t629 = 0.2e1 * t675 - t805 * m(2) - Icges(2,3) - Icges(3,3) + (-0.2e1 * t726 * rSges(3,1) * pkin(1) - t759 - t761 - t764) * m(3);
t771 = t648 * ((t599 * t857 * t848 + (t663 * t841 - t821 + t827 / 0.2e1 + t818 + t780 * t846) * t596) * t596 + t590 * t602 - t629 * t796);
t770 = t649 * ((t600 * t855 * t848 + (t664 * t841 - t820 + t826 / 0.2e1 + t817 + t779 * t846) * t597) * t597 + t591 * t603 - t629 * t795);
t769 = t650 * ((t601 * t853 * t848 + (t665 * t841 - t819 + t825 / 0.2e1 + t816 + t778 * t846) * t598) * t598 + t592 * t604 - t629 * t794);
t1 = [t697 * t769 + t696 * t770 + t695 * t771 + (t788 * t615 - t622 * t791) * t718 + (t789 * t613 - t621 * t792) * t717 + (t790 * t611 - t620 * t793) * t716; t700 * t769 + t699 * t770 + t698 * t771 + (t788 * t616 - t619 * t791) * t718 + (t789 * t614 - t618 * t792) * t717 + (t790 * t612 - t617 * t793) * t716; (t788 * t638 + t809 * t747) * t718 + (t789 * t637 + t810 * t745) * t717 + (t790 * t636 + t811 * t743) * t716;];
taucX  = t1;
