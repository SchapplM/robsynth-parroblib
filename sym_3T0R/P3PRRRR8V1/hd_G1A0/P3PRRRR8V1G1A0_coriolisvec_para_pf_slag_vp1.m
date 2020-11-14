% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 16:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:49:44
% EndTime: 2020-08-06 16:49:48
% DurationCPUTime: 3.54s
% Computational Cost: add. (13248->248), mult. (35478->507), div. (3729->9), fcn. (42552->28), ass. (0->229)
t741 = cos(qJ(3,2));
t723 = 0.1e1 / t741;
t742 = cos(qJ(2,2));
t736 = sin(qJ(2,2));
t798 = t736 * t741;
t696 = pkin(2) * t798 - t742 * pkin(5);
t727 = sin(pkin(3));
t729 = cos(pkin(3));
t735 = sin(qJ(3,2));
t840 = pkin(2) * t735;
t771 = 0.1e1 / (t696 * t727 + t729 * t840);
t820 = t771 * t723;
t737 = sin(qJ(3,1));
t743 = cos(qJ(3,1));
t768 = rSges(3,1) * t743 - rSges(3,2) * t737;
t860 = t768 * m(3);
t769 = rSges(3,1) * t741 - rSges(3,2) * t735;
t859 = t769 * m(3);
t733 = sin(qJ(3,3));
t739 = cos(qJ(3,3));
t770 = rSges(3,1) * t739 - rSges(3,2) * t733;
t858 = t770 * m(3);
t730 = legFrame(3,3);
t713 = sin(t730);
t716 = cos(t730);
t726 = sin(pkin(6));
t728 = cos(pkin(6));
t683 = -t726 * t713 + t716 * t728;
t740 = cos(qJ(2,3));
t734 = sin(qJ(2,3));
t808 = t729 * t734;
t689 = t726 * t740 + t728 * t808;
t692 = -t726 * t808 + t728 * t740;
t815 = t727 * t739;
t649 = (-t689 * t716 - t713 * t692) * t733 - t683 * t815;
t686 = t728 * t713 + t716 * t726;
t652 = (-t713 * t689 + t692 * t716) * t733 - t686 * t815;
t801 = t734 * t739;
t695 = pkin(2) * t801 - t740 * pkin(5);
t809 = t729 * t733;
t677 = 0.1e1 / (pkin(2) * t809 + t695 * t727);
t721 = 0.1e1 / t739;
t748 = xDP(2);
t749 = xDP(1);
t634 = (t649 * t749 + t652 * t748) * t721 * t677;
t720 = t739 ^ 2;
t755 = pkin(5) ^ 2;
t756 = pkin(2) ^ 2;
t793 = t740 * t686;
t794 = t740 * t683;
t802 = t734 * t686;
t803 = t734 * t683;
t839 = pkin(2) * t739;
t643 = -(t729 * t794 - t802) * t839 - pkin(5) * (t729 * t803 + t793);
t644 = -(t729 * t793 + t803) * t839 - (t729 * t802 - t794) * pkin(5);
t757 = 0.1e1 / pkin(2);
t774 = t727 * t801;
t814 = t727 * t740;
t854 = 0.1e1 / (-pkin(5) * t814 + (t774 + t809) * pkin(2));
t821 = t854 * t721;
t628 = (t643 * t749 + t644 * t748) * t757 * t821;
t846 = pkin(2) * t628;
t783 = t733 * t846;
t616 = -pkin(5) * t783 + (t720 * t756 + t755) * t634;
t836 = pkin(5) * t634;
t780 = t733 * t836;
t619 = t780 - t846;
t610 = (t616 * t634 - t619 * t846) * t854;
t625 = t628 ^ 2;
t631 = t634 ^ 2;
t701 = t733 * rSges(3,1) + t739 * rSges(3,2);
t851 = m(3) * rSges(3,3);
t709 = m(2) * rSges(2,2) - t851;
t719 = -m(1) - m(2) - m(3);
t747 = m(2) * rSges(2,1);
t784 = 0.2e1 * m(3);
t857 = ((-t631 * t747 - (t631 + t625) * t858) * t734 - (t701 * t628 * t784 + t634 * t709) * t740 * t634) * t727 - t719 * t610;
t731 = legFrame(2,3);
t714 = sin(t731);
t717 = cos(t731);
t684 = -t726 * t714 + t717 * t728;
t807 = t729 * t736;
t690 = t726 * t742 + t728 * t807;
t693 = -t726 * t807 + t728 * t742;
t813 = t727 * t741;
t650 = (-t690 * t717 - t714 * t693) * t735 - t684 * t813;
t687 = t728 * t714 + t717 * t726;
t653 = (-t714 * t690 + t693 * t717) * t735 - t687 * t813;
t635 = (t650 * t749 + t653 * t748) * t820;
t835 = pkin(5) * t635;
t779 = t735 * t835;
t791 = t742 * t687;
t792 = t742 * t684;
t799 = t736 * t687;
t800 = t736 * t684;
t838 = pkin(2) * t741;
t645 = -(t729 * t792 - t799) * t838 - pkin(5) * (t729 * t800 + t791);
t646 = -(t729 * t791 + t800) * t838 - (t729 * t799 - t792) * pkin(5);
t629 = (t645 * t749 + t646 * t748) * t757 * t820;
t845 = pkin(2) * t629;
t620 = t779 - t845;
t722 = t741 ^ 2;
t782 = t629 * t840;
t830 = (-pkin(5) * t782 + (t722 * t756 + t755) * t635) * t635;
t611 = (t620 * t845 - t830) * t771;
t626 = t629 ^ 2;
t632 = t635 ^ 2;
t702 = t735 * rSges(3,1) + t741 * rSges(3,2);
t856 = ((-t632 * t747 - (t632 + t626) * t859) * t736 - (t702 * t629 * t784 + t635 * t709) * t742 * t635) * t727 + t719 * t611;
t732 = legFrame(1,3);
t715 = sin(t732);
t718 = cos(t732);
t685 = -t726 * t715 + t718 * t728;
t744 = cos(qJ(2,1));
t738 = sin(qJ(2,1));
t805 = t729 * t738;
t691 = t726 * t744 + t728 * t805;
t694 = -t726 * t805 + t728 * t744;
t811 = t727 * t743;
t651 = (-t691 * t718 - t715 * t694) * t737 - t685 * t811;
t688 = t728 * t715 + t718 * t726;
t654 = (-t715 * t691 + t694 * t718) * t737 - t688 * t811;
t795 = t738 * t743;
t697 = pkin(2) * t795 - t744 * pkin(5);
t806 = t729 * t737;
t679 = 0.1e1 / (pkin(2) * t806 + t697 * t727);
t725 = 0.1e1 / t743;
t636 = (t651 * t749 + t654 * t748) * t725 * t679;
t724 = t743 ^ 2;
t789 = t744 * t688;
t790 = t744 * t685;
t796 = t738 * t688;
t797 = t738 * t685;
t837 = pkin(2) * t743;
t647 = -(t729 * t790 - t796) * t837 - pkin(5) * (t729 * t797 + t789);
t648 = -(t729 * t789 + t797) * t837 - (t729 * t796 - t790) * pkin(5);
t772 = t727 * t795;
t810 = t727 * t744;
t853 = 0.1e1 / (-pkin(5) * t810 + (t772 + t806) * pkin(2));
t819 = t853 * t725;
t630 = (t647 * t749 + t648 * t748) * t757 * t819;
t844 = pkin(2) * t630;
t781 = t737 * t844;
t618 = -pkin(5) * t781 + (t724 * t756 + t755) * t636;
t834 = pkin(5) * t636;
t778 = t737 * t834;
t621 = -t778 + t844;
t612 = (t618 * t636 + t621 * t844) * t853;
t627 = t630 ^ 2;
t633 = t636 ^ 2;
t703 = t737 * rSges(3,1) + t743 * rSges(3,2);
t855 = ((-t633 * t747 - (t633 + t627) * t860) * t738 - (t703 * t630 * t784 + t636 * t709) * t744 * t636) * t727 - t719 * t612;
t712 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t852 = 0.2e1 * t712;
t753 = rSges(3,2) ^ 2;
t754 = rSges(3,1) ^ 2;
t705 = (-t753 + t754) * m(3) - Icges(3,1) + Icges(3,2);
t850 = -t705 / 0.2e1;
t710 = rSges(3,2) * t851 - Icges(3,6);
t849 = -t710 / 0.4e1;
t711 = rSges(3,1) * t851 - Icges(3,5);
t848 = t711 / 0.4e1;
t847 = m(3) * t729;
t843 = pkin(2) * t720;
t842 = pkin(2) * t722;
t841 = pkin(2) * t724;
t804 = t729 * t757;
t818 = t727 * t733;
t826 = t634 * t854;
t604 = (t616 * t804 * t826 + (-t628 * t695 * t818 + t729 * (t628 * t843 - t780)) * t677 * t628) * t721;
t664 = -t727 * t734 * t701 + t770 * t729;
t833 = t604 * t664;
t817 = t727 * t735;
t605 = (t804 * t830 + (-t629 * t696 * t817 + t729 * (t629 * t842 - t779)) * t629) * t820;
t665 = -t727 * t736 * t702 + t769 * t729;
t832 = t605 * t665;
t816 = t727 * t737;
t825 = t636 * t853;
t606 = (t618 * t804 * t825 + (-t630 * t697 * t816 + t729 * (t630 * t841 - t778)) * t679 * t630) * t725;
t666 = -t727 * t738 * t703 + t768 * t729;
t831 = t606 * t666;
t829 = t625 * t701;
t828 = t626 * t702;
t827 = t627 * t703;
t667 = (t747 + t858) * t740 - t734 * t709;
t824 = t667 * t727;
t668 = (t747 + t859) * t742 - t736 * t709;
t823 = t668 * t727;
t669 = (t747 + t860) * t744 - t738 * t709;
t822 = t669 * t727;
t812 = t727 * t742;
t595 = (((t729 * t628 + t634 * t814) * t843 - (t783 - t836) * t774 + t729 * t619) * t826 + (t628 * t814 + (t720 * t729 - t733 * t774 - t729) * t634) * t854 * t846) * t721;
t788 = -m(3) * t833 - t595 * t824 - t829 * t847 + t857;
t773 = t727 * t798;
t596 = (((t729 * t629 + t635 * t812) * t842 - (t782 - t835) * t773 + t729 * t620) * t635 + (t629 * t812 + (t722 * t729 - t735 * t773 - t729) * t635) * t845) * t820;
t787 = -m(3) * t832 - t596 * t823 - t828 * t847 + t856;
t597 = (((t729 * t630 + t636 * t810) * t841 - (t781 - t834) * t772 - t729 * t621) * t825 + (t630 * t810 + (t724 * t729 - t737 * t772 - t729) * t636) * t853 * t844) * t725;
t786 = -m(3) * t831 - t597 * t822 - t827 * t847 + t855;
t785 = t753 + t754;
t777 = t705 * t733 * t739;
t776 = t705 * t735 * t741;
t775 = t705 * t737 * t743;
t767 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - (0.2e1 * rSges(3,3) ^ 2 + t785) * m(3) / 0.2e1 - Icges(3,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,3);
t680 = t710 * t739 + t711 * t733;
t750 = 0.2e1 * qJ(3,3);
t766 = (-0.4e1 * ((t733 * t849 + t739 * t848) * t628 + (t777 / 0.2e1 + (t720 - 0.1e1 / 0.2e1) * t712) * t634) * t628 + t610 * t824 + (cos(t750) * t850 + t712 * sin(t750) + t767) * t595 + t680 * t604) * t721;
t681 = t710 * t741 + t711 * t735;
t751 = 0.2e1 * qJ(3,2);
t765 = (-0.4e1 * ((t735 * t849 + t741 * t848) * t629 + (t776 / 0.2e1 + (t722 - 0.1e1 / 0.2e1) * t712) * t635) * t629 - t611 * t823 + (cos(t751) * t850 + t712 * sin(t751) + t767) * t596 + t681 * t605) * t723;
t682 = t710 * t743 + t711 * t737;
t752 = 0.2e1 * qJ(3,1);
t764 = (-0.4e1 * ((t737 * t849 + t743 * t848) * t630 + (t775 / 0.2e1 + (t724 - 0.1e1 / 0.2e1) * t712) * t636) * t630 + t612 * t822 + (cos(t752) * t850 + t712 * sin(t752) + t767) * t597 + t682 * t606) * t725;
t763 = pkin(2) * t818 - t695 * t729;
t762 = pkin(2) * t817 - t696 * t729;
t761 = pkin(2) * t816 - t697 * t729;
t707 = -t785 * m(3) - Icges(3,3);
t760 = (t664 * m(3) * t610 + t680 * t595 + t707 * t604 + t631 * (t720 * t852 - t712 + t777)) * t821;
t759 = (-t665 * m(3) * t611 + t681 * t596 + t707 * t605 + t632 * (t722 * t852 - t712 + t776)) * t820;
t758 = (t666 * m(3) * t612 + t682 * t597 + t707 * t606 + t633 * (t724 * t852 - t712 + t775)) * t819;
t700 = pkin(5) * t738 + t744 * t837;
t699 = pkin(5) * t736 + t742 * t838;
t698 = pkin(5) * t734 + t740 * t839;
t660 = t726 * t700 - t761 * t728;
t659 = t726 * t699 - t762 * t728;
t658 = t726 * t698 - t763 * t728;
t657 = t700 * t728 + t761 * t726;
t656 = t699 * t728 + t762 * t726;
t655 = t698 * t728 + t763 * t726;
t1 = [(t651 * t764 + t786 * (t657 * t718 - t715 * t660)) * t679 + (t650 * t765 + t787 * (t656 * t717 - t714 * t659)) * t771 + (t649 * t766 + t788 * (t655 * t716 - t713 * t658)) * t677 + (t643 * t760 + t645 * t759 + t647 * t758) * t757; (t654 * t764 + t786 * (t657 * t715 + t660 * t718)) * t679 + (t653 * t765 + t787 * (t656 * t714 + t659 * t717)) * t771 + (t652 * t766 + t788 * (t655 * t713 + t658 * t716)) * t677 + (t644 * t760 + t646 * t759 + t648 * t758) * t757; (-t595 * t667 - t596 * t668 - t597 * t669) * t727 + (-t833 - t832 - t831 + (-t827 - t828 - t829) * t729) * m(3) + t855 + t856 + t857;];
taucX  = t1;
