% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR9V1G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:48
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:46:50
% EndTime: 2020-08-06 18:46:54
% DurationCPUTime: 3.66s
% Computational Cost: add. (16113->334), mult. (14238->548), div. (3540->11), fcn. (11907->55), ass. (0->258)
t703 = (pkin(7) + qJ(3,1));
t656 = 2 * t703;
t650 = cos(t656);
t884 = (t650 + 0.1e1) * pkin(3);
t702 = (pkin(7) + qJ(3,2));
t655 = 2 * t702;
t649 = cos(t655);
t883 = (t649 + 0.1e1) * pkin(3);
t701 = (pkin(7) + qJ(3,3));
t654 = 2 * t701;
t648 = cos(t654);
t882 = (t648 + 0.1e1) * pkin(3);
t881 = 2 * pkin(1);
t880 = -pkin(3) / 0.2e1;
t736 = 2 * pkin(7);
t700 = t736 + qJ(3,1);
t673 = cos(t700);
t720 = cos(qJ(3,1));
t810 = t673 + t720;
t873 = t810 * pkin(2);
t879 = t884 / 0.2e1 + t873 / 0.2e1;
t699 = t736 + qJ(3,2);
t672 = cos(t699);
t718 = cos(qJ(3,2));
t811 = t672 + t718;
t872 = t811 * pkin(2);
t878 = t883 / 0.2e1 + t872 / 0.2e1;
t698 = t736 + qJ(3,3);
t671 = cos(t698);
t716 = cos(qJ(3,3));
t812 = t671 + t716;
t871 = t812 * pkin(2);
t877 = t882 / 0.2e1 + t871 / 0.2e1;
t726 = xDP(3);
t704 = t726 ^ 2;
t745 = 0.1e1 / pkin(3);
t870 = t704 * t745;
t846 = pkin(5) + qJ(2,3);
t692 = -pkin(6) - t846;
t684 = 0.1e1 / t692;
t645 = sin(t654);
t668 = sin(t701);
t665 = sin(t698);
t710 = sin(qJ(3,3));
t815 = t665 + t710;
t614 = pkin(2) * t815 + pkin(3) * t645 + t668 * t881;
t674 = cos(t701);
t659 = 0.1e1 / t674;
t792 = t614 * t659 / 0.2e1;
t706 = cos(pkin(7));
t677 = t706 * pkin(2);
t653 = t677 + pkin(1);
t711 = sin(qJ(1,3));
t717 = cos(qJ(1,3));
t617 = t653 * t711 + t692 * t717;
t620 = t653 * t717 - t692 * t711;
t707 = legFrame(3,3);
t678 = sin(t707);
t681 = cos(t707);
t627 = t678 * t717 + t681 * t711;
t854 = pkin(3) * t674;
t599 = t617 * t681 + t620 * t678 + t627 * t854;
t727 = xDP(2);
t842 = t599 * t727;
t624 = -t678 * t711 + t681 * t717;
t596 = -t617 * t678 + t620 * t681 + t624 * t854;
t728 = xDP(1);
t845 = t596 * t728;
t590 = (t726 * t792 + t842 + t845) * t684;
t831 = t659 * t726;
t605 = (t624 * t728 + t627 * t727 + t668 * t831) * t684;
t602 = pkin(1) * t605;
t735 = 3 * pkin(7);
t737 = 2 * qJ(3,3);
t744 = pkin(3) ^ 2;
t747 = pkin(2) ^ 2;
t821 = t745 * t726;
t797 = t659 * t821;
t780 = t645 * t797;
t748 = pkin(1) ^ 2;
t793 = -0.3e1 * t744 - 0.2e1 * t747 - (4 * t748);
t794 = -0.4e1 * pkin(3) * t677;
t803 = t614 * t831;
t819 = -0.2e1 * pkin(2) * pkin(3);
t832 = t659 * t684;
t838 = t605 * t684;
t855 = -t726 / 0.2e1;
t866 = 0.1e1 / t674 ^ 2;
t575 = (t653 + t854) * t866 * t870 * t832 + (t692 * t645 * t866 * t855 + (-(-t744 * cos((3 * t701)) + (-0.4e1 * t692 ^ 2 + t793) * t674 + t794 + (cos((t737 + t735)) + cos((pkin(7) + t737))) * t819 + (-cos((-pkin(7) + qJ(3,3))) - cos((t735 + qJ(3,3)))) * t747) * t605 / 0.4e1 + t692 * t780 * t880 + (-t871 - t882) * (-t602 - (-t845 / 0.2e1 - t842 / 0.2e1 - t803 / 0.4e1) * t684) - (t674 * t881 + t877) * t590) * t659) * t838;
t660 = t659 * t866;
t758 = t803 + 0.2e1 * t842 + 0.2e1 * t845;
t789 = -t832 / 0.2e1;
t839 = t605 * t659;
t581 = (-t660 * t870 - (-t590 + t839 * t877 + (t602 * t659 + t758 * t789) * t674) * t605) * t684;
t729 = (m(2) + m(3));
t857 = rSges(2,2) * m(2);
t807 = sin(pkin(7)) * t857;
t861 = pkin(2) * m(3);
t833 = (m(2) * rSges(2,1) + t861) * t706;
t762 = (pkin(1) * t729) - t807 + t833;
t769 = rSges(3,1) * t674 - rSges(3,2) * t668;
t608 = t769 * m(3) + t762;
t689 = rSges(3,3) + t846;
t851 = t689 * m(3);
t638 = -rSges(3,2) * t851 + Icges(3,6);
t641 = rSges(3,1) * t851 - Icges(3,5);
t611 = -t638 * t674 + t641 * t668;
t632 = t851 + m(2) * (rSges(2,3) + qJ(2,3));
t740 = rSges(3,2) ^ 2;
t742 = rSges(3,1) ^ 2;
t636 = m(3) * (-t740 + t742) - Icges(3,1) + Icges(3,2);
t858 = m(3) * rSges(3,2);
t658 = rSges(3,1) * t858 - Icges(3,4);
t741 = rSges(2,2) ^ 2;
t743 = rSges(2,1) ^ 2;
t869 = -2 * pkin(1);
t755 = -(m(3) * t747 + (-t741 + t743) * m(2) - Icges(2,1) + Icges(2,2)) * cos(t736) / 0.2e1 - (-rSges(2,1) * t857 + Icges(2,4)) * sin(t736) + t833 * t869 - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t807 * t881 - Icges(3,2) / 0.2e1 - Icges(2,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,1) / 0.2e1 - Icges(1,3);
t733 = 2 * t748;
t763 = -t733 / 0.2e1 - t740 / 0.2e1 - t742 / 0.2e1 - t747 / 0.2e1;
t806 = m(3) * t821;
t785 = pkin(2) * t806;
t772 = t659 * t785;
t774 = 0.2e1 * rSges(2,3) ^ 2 + t733 + t741 + t743;
t825 = t704 / pkin(3) ^ 2;
t777 = t660 * t668 * t825;
t784 = -0.2e1 * t658 * t821;
t804 = t858 * t869;
t805 = m(3) * rSges(3,1) * t881;
t856 = -t636 / 0.2e1;
t859 = -rSges(3,2) / 0.2e1;
t860 = -rSges(3,1) / 0.2e1;
t862 = -m(2) / 0.2e1;
t867 = 0.4e1 * rSges(2,3);
t818 = t608 * t575 - t611 * t777 - (t648 * t856 + t658 * t645 + ((t867 + 0.2e1 * qJ(2,3)) * qJ(2,3) + t774) * t862 + t755 + (-t689 ^ 2 - t769 * t881 + (-rSges(3,1) * t812 + rSges(3,2) * t815) * pkin(2) + t763) * m(3)) * t581 + t605 * t636 * t780 - t648 * t784 * t839 + (-t605 * t804 - t641 * t797) * t821 - (-t605 * t805 + t638 * t797) * t668 * t797 - 0.2e1 * t605 * (-t590 * t632 + (t710 * t860 + t716 * t859) * t772) - (-rSges(3,1) * t665 - rSges(3,2) * t671) * t605 * t772;
t847 = pkin(5) + qJ(2,2);
t693 = -pkin(6) - t847;
t685 = 0.1e1 / t693;
t646 = sin(t655);
t669 = sin(t702);
t666 = sin(t699);
t712 = sin(qJ(3,2));
t814 = t666 + t712;
t615 = pkin(2) * t814 + pkin(3) * t646 + t669 * t881;
t675 = cos(t702);
t661 = 0.1e1 / t675;
t791 = t615 * t661 / 0.2e1;
t713 = sin(qJ(1,2));
t719 = cos(qJ(1,2));
t618 = t653 * t713 + t693 * t719;
t621 = t653 * t719 - t693 * t713;
t708 = legFrame(2,3);
t679 = sin(t708);
t682 = cos(t708);
t628 = t679 * t719 + t682 * t713;
t853 = pkin(3) * t675;
t600 = t618 * t682 + t621 * t679 + t628 * t853;
t841 = t600 * t727;
t625 = -t679 * t713 + t682 * t719;
t597 = -t618 * t679 + t621 * t682 + t625 * t853;
t844 = t597 * t728;
t591 = (t726 * t791 + t841 + t844) * t685;
t829 = t661 * t726;
t606 = (t625 * t728 + t628 * t727 + t669 * t829) * t685;
t603 = pkin(1) * t606;
t738 = 2 * qJ(3,2);
t796 = t661 * t821;
t779 = t646 * t796;
t802 = t615 * t829;
t830 = t661 * t685;
t836 = t606 * t685;
t865 = 0.1e1 / t675 ^ 2;
t576 = (t653 + t853) * t865 * t870 * t830 + (t693 * t646 * t865 * t855 + (-(-t744 * cos((3 * t702)) + (-0.4e1 * t693 ^ 2 + t793) * t675 + t794 + (cos((t735 + t738)) + cos((pkin(7) + t738))) * t819 + (-cos((-pkin(7) + qJ(3,2))) - cos((t735 + qJ(3,2)))) * t747) * t606 / 0.4e1 + t693 * t779 * t880 + (-t872 - t883) * (-t603 - (-t844 / 0.2e1 - t841 / 0.2e1 - t802 / 0.4e1) * t685) - (t675 * t881 + t878) * t591) * t661) * t836;
t662 = t661 * t865;
t757 = t802 + 0.2e1 * t841 + 0.2e1 * t844;
t788 = -t830 / 0.2e1;
t837 = t606 * t661;
t582 = (-t662 * t870 - (-t591 + t837 * t878 + (t603 * t661 + t757 * t788) * t675) * t606) * t685;
t768 = rSges(3,1) * t675 - rSges(3,2) * t669;
t609 = t768 * m(3) + t762;
t690 = rSges(3,3) + t847;
t850 = t690 * m(3);
t639 = -rSges(3,2) * t850 + Icges(3,6);
t642 = rSges(3,1) * t850 - Icges(3,5);
t612 = -t639 * t675 + t642 * t669;
t633 = t850 + m(2) * (rSges(2,3) + qJ(2,2));
t771 = t661 * t785;
t776 = t662 * t669 * t825;
t817 = t609 * t576 - t612 * t776 - (t649 * t856 + t658 * t646 + ((t867 + 0.2e1 * qJ(2,2)) * qJ(2,2) + t774) * t862 + t755 + (-t690 ^ 2 - t768 * t881 + (-rSges(3,1) * t811 + rSges(3,2) * t814) * pkin(2) + t763) * m(3)) * t582 + t606 * t636 * t779 - t649 * t784 * t837 + (-t606 * t804 - t642 * t796) * t821 - (-t606 * t805 + t639 * t796) * t669 * t796 - 0.2e1 * t606 * (-t591 * t633 + (t712 * t860 + t718 * t859) * t771) - (-rSges(3,1) * t666 - rSges(3,2) * t672) * t606 * t771;
t848 = pkin(5) + qJ(2,1);
t694 = -pkin(6) - t848;
t686 = 0.1e1 / t694;
t647 = sin(t656);
t670 = sin(t703);
t667 = sin(t700);
t714 = sin(qJ(3,1));
t813 = t667 + t714;
t616 = pkin(2) * t813 + pkin(3) * t647 + t670 * t881;
t676 = cos(t703);
t663 = 0.1e1 / t676;
t790 = t616 * t663 / 0.2e1;
t715 = sin(qJ(1,1));
t721 = cos(qJ(1,1));
t619 = t653 * t715 + t694 * t721;
t622 = t653 * t721 - t694 * t715;
t709 = legFrame(1,3);
t680 = sin(t709);
t683 = cos(t709);
t629 = t680 * t721 + t683 * t715;
t852 = pkin(3) * t676;
t601 = t619 * t683 + t622 * t680 + t629 * t852;
t840 = t601 * t727;
t626 = -t680 * t715 + t683 * t721;
t598 = -t619 * t680 + t622 * t683 + t626 * t852;
t843 = t598 * t728;
t592 = (t726 * t790 + t840 + t843) * t686;
t827 = t663 * t726;
t607 = (t626 * t728 + t629 * t727 + t670 * t827) * t686;
t604 = pkin(1) * t607;
t739 = 2 * qJ(3,1);
t795 = t663 * t821;
t778 = t647 * t795;
t801 = t616 * t827;
t828 = t663 * t686;
t834 = t607 * t686;
t864 = 0.1e1 / t676 ^ 2;
t577 = (t653 + t852) * t864 * t870 * t828 + (t694 * t647 * t864 * t855 + (-(-t744 * cos((3 * t703)) + (-0.4e1 * t694 ^ 2 + t793) * t676 + t794 + (cos((t735 + t739)) + cos((pkin(7) + t739))) * t819 + (-cos((-pkin(7) + qJ(3,1))) - cos((t735 + qJ(3,1)))) * t747) * t607 / 0.4e1 + t694 * t778 * t880 + (-t873 - t884) * (-t604 - (-t843 / 0.2e1 - t840 / 0.2e1 - t801 / 0.4e1) * t686) - (t676 * t881 + t879) * t592) * t663) * t834;
t664 = t663 * t864;
t756 = t801 + 0.2e1 * t840 + 0.2e1 * t843;
t787 = -t828 / 0.2e1;
t835 = t607 * t663;
t583 = (-t664 * t870 - (-t592 + t835 * t879 + (t604 * t663 + t756 * t787) * t676) * t607) * t686;
t767 = rSges(3,1) * t676 - rSges(3,2) * t670;
t610 = t767 * m(3) + t762;
t691 = rSges(3,3) + t848;
t849 = t691 * m(3);
t640 = -rSges(3,2) * t849 + Icges(3,6);
t643 = rSges(3,1) * t849 - Icges(3,5);
t613 = -t640 * t676 + t643 * t670;
t634 = t849 + m(2) * (rSges(2,3) + qJ(2,1));
t770 = t663 * t785;
t775 = t664 * t670 * t825;
t816 = t610 * t577 - t613 * t775 - (t650 * t856 + t658 * t647 + ((t867 + 0.2e1 * qJ(2,1)) * qJ(2,1) + t774) * t862 + t755 + (-t691 ^ 2 - t767 * t881 + (-rSges(3,1) * t810 + rSges(3,2) * t813) * pkin(2) + t763) * m(3)) * t583 + t607 * t636 * t778 - t650 * t784 * t835 + (-t607 * t804 - t643 * t795) * t821 - (-t607 * t805 + t640 * t795) * t670 * t795 - 0.2e1 * t607 * (-t592 * t634 + (t714 * t860 + t720 * t859) * t770) - (-rSges(3,1) * t667 - rSges(3,2) * t673) * t607 * t770;
t863 = 0.2e1 * t658;
t809 = 0.2e1 * m(3);
t786 = -0.2e1 * t806;
t572 = -t575 * t729 - t581 * t608;
t652 = rSges(3,2) * t786;
t773 = rSges(3,1) * t786;
t593 = t659 * t668 * t773 - t605 * t632 + t652;
t783 = t605 * t593 + t572;
t573 = -t576 * t729 - t582 * t609;
t594 = t661 * t669 * t773 - t606 * t633 + t652;
t782 = t606 * t594 + t573;
t574 = -t577 * t729 - t583 * t610;
t595 = t663 * t670 * t773 - t607 * t634 + t652;
t781 = t607 * t595 + t574;
t644 = -(t740 + t742) * m(3) - Icges(3,3);
t1 = [-(t598 * t781 + t626 * t816) * t686 - (t597 * t782 + t625 * t817) * t685 - (t596 * t783 + t624 * t818) * t684; -(t601 * t781 + t629 * t816) * t686 - (t600 * t782 + t628 * t817) * t685 - (t599 * t783 + t627 * t818) * t684; t614 * t572 * t789 + t615 * t573 * t788 + t616 * t574 * t787 - t593 * t792 * t838 - t594 * t791 * t836 - t595 * t790 * t834 - t818 * t668 * t832 - t817 * t669 * t830 - t816 * t670 * t828 + ((-t581 * t611 - t644 * t777 - t605 * ((rSges(3,1) * t668 + rSges(3,2) * t674) * (t758 * t684 - t602) * t809 - (t636 * t645 + t648 * t863 + (rSges(3,1) * t815 + rSges(3,2) * t812) * t861) * t605) / 0.2e1) * t659 + (-t582 * t612 - t644 * t776 - t606 * ((rSges(3,1) * t669 + rSges(3,2) * t675) * (t685 * t757 - t603) * t809 - (t636 * t646 + t649 * t863 + (rSges(3,1) * t814 + rSges(3,2) * t811) * t861) * t606) / 0.2e1) * t661 + (-t583 * t613 - t644 * t775 - t607 * ((rSges(3,1) * t670 + rSges(3,2) * t676) * (t686 * t756 - t604) * t809 - (t636 * t647 + t650 * t863 + (rSges(3,1) * t813 + rSges(3,2) * t810) * t861) * t607) / 0.2e1) * t663) * t745;];
taucX  = t1;
