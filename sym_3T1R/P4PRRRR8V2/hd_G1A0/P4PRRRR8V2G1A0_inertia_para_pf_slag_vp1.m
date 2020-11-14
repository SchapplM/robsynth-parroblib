% Calculate inertia matrix for parallel robot
% P4PRRRR8V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% MX [4x4]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRRR8V2G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:09:32
% EndTime: 2020-08-07 11:09:35
% DurationCPUTime: 3.20s
% Computational Cost: add. (13609->387), mult. (26298->701), div. (896->9), fcn. (25640->30), ass. (0->299)
t786 = cos(qJ(2,1));
t780 = sin(qJ(2,1));
t789 = pkin(7) + pkin(6);
t836 = t780 * t789;
t730 = pkin(2) * t786 + t836;
t763 = sin(pkin(8));
t765 = cos(pkin(8));
t745 = t789 * t786;
t727 = pkin(2) * t780 - t745;
t766 = cos(pkin(4));
t764 = sin(pkin(4));
t779 = sin(qJ(3,1));
t852 = t764 * t779;
t805 = pkin(3) * t852 - t727 * t766;
t903 = t730 * t765 + t763 * t805;
t784 = cos(qJ(2,2));
t778 = sin(qJ(2,2));
t837 = t778 * t789;
t729 = pkin(2) * t784 + t837;
t744 = t789 * t784;
t726 = pkin(2) * t778 - t744;
t777 = sin(qJ(3,2));
t854 = t764 * t777;
t806 = pkin(3) * t854 - t726 * t766;
t902 = t729 * t765 + t763 * t806;
t782 = cos(qJ(2,3));
t776 = sin(qJ(2,3));
t838 = t776 * t789;
t728 = pkin(2) * t782 + t838;
t743 = t789 * t782;
t725 = pkin(2) * t776 - t743;
t775 = sin(qJ(3,3));
t856 = t764 * t775;
t807 = pkin(3) * t856 - t725 * t766;
t901 = t728 * t765 + t763 * t807;
t774 = cos(qJ(2,4));
t772 = sin(qJ(2,4));
t839 = t772 * t789;
t724 = pkin(2) * t774 + t839;
t738 = t789 * t774;
t723 = pkin(2) * t772 - t738;
t771 = sin(qJ(3,4));
t859 = t764 * t771;
t808 = pkin(3) * t859 - t723 * t766;
t900 = t724 * t765 + t763 * t808;
t898 = (m(3) * rSges(3,2));
t899 = -2 * rSges(3,1) * t898 + 2 * Icges(3,4);
t790 = pkin(2) * m(3);
t773 = cos(qJ(3,4));
t812 = rSges(3,1) * t773 - rSges(3,2) * t771;
t858 = t764 * t772;
t897 = m(3) * (t812 * t766 - (t771 * rSges(3,1) + t773 * rSges(3,2)) * t858);
t781 = cos(qJ(3,3));
t811 = rSges(3,1) * t781 - rSges(3,2) * t775;
t855 = t764 * t776;
t896 = m(3) * (t811 * t766 - (t775 * rSges(3,1) + t781 * rSges(3,2)) * t855);
t783 = cos(qJ(3,2));
t810 = rSges(3,1) * t783 - rSges(3,2) * t777;
t853 = t764 * t778;
t895 = m(3) * (t810 * t766 - (t777 * rSges(3,1) + t783 * rSges(3,2)) * t853);
t785 = cos(qJ(3,1));
t809 = rSges(3,1) * t785 - rSges(3,2) * t779;
t851 = t764 * t780;
t894 = m(3) * (t809 * t766 - (t779 * rSges(3,1) + t785 * rSges(3,2)) * t851);
t758 = t773 ^ 2;
t893 = pkin(3) * t758;
t760 = t781 ^ 2;
t892 = pkin(3) * t760;
t761 = t783 ^ 2;
t891 = pkin(3) * t761;
t762 = t785 ^ 2;
t890 = pkin(3) * t762;
t787 = pkin(6) + rSges(3,3);
t889 = t787 * m(3);
t767 = legFrame(4,3);
t747 = sin(t767);
t751 = cos(t767);
t689 = -t763 * t747 + t751 * t765;
t693 = t765 * t747 + t751 * t763;
t737 = t773 * pkin(3) + pkin(2);
t706 = t772 * t737 - t738;
t868 = (t737 * t774 + t839) * t766;
t653 = -t689 * t868 + t706 * t693;
t847 = t766 * t771;
t857 = t764 * t773;
t673 = 0.1e1 / (t706 * t857 + t737 * t847);
t888 = t653 * t673;
t654 = -t706 * t689 - t693 * t868;
t887 = t654 * t673;
t768 = legFrame(3,3);
t748 = sin(t768);
t752 = cos(t768);
t690 = -t763 * t748 + t752 * t765;
t694 = t765 * t748 + t752 * t763;
t739 = t781 * pkin(3) + pkin(2);
t712 = t776 * t739 - t743;
t867 = (t739 * t782 + t838) * t766;
t655 = -t690 * t867 + t712 * t694;
t845 = t766 * t775;
t850 = t764 * t781;
t674 = 0.1e1 / (t712 * t850 + t739 * t845);
t886 = t655 * t674;
t769 = legFrame(2,3);
t749 = sin(t769);
t753 = cos(t769);
t691 = -t763 * t749 + t753 * t765;
t695 = t765 * t749 + t753 * t763;
t740 = t783 * pkin(3) + pkin(2);
t713 = t778 * t740 - t744;
t866 = (t740 * t784 + t837) * t766;
t656 = -t691 * t866 + t713 * t695;
t843 = t766 * t777;
t849 = t764 * t783;
t675 = 0.1e1 / (t713 * t849 + t740 * t843);
t885 = t656 * t675;
t770 = legFrame(1,3);
t750 = sin(t770);
t754 = cos(t770);
t692 = -t763 * t750 + t754 * t765;
t696 = t765 * t750 + t754 * t763;
t741 = t785 * pkin(3) + pkin(2);
t714 = t780 * t741 - t745;
t865 = (t741 * t786 + t836) * t766;
t657 = -t692 * t865 + t714 * t696;
t841 = t766 * t779;
t848 = t764 * t785;
t676 = 0.1e1 / (t714 * t848 + t741 * t841);
t884 = t657 * t676;
t658 = -t712 * t690 - t694 * t867;
t883 = t658 * t674;
t659 = -t713 * t691 - t695 * t866;
t882 = t659 * t675;
t660 = -t714 * t692 - t696 * t865;
t881 = t660 * t676;
t665 = 0.1e1 / (t858 * t893 + (pkin(3) * t847 + t723 * t764) * t773 + pkin(2) * t847);
t759 = m(1) + m(2) + m(3);
t880 = t665 * t759;
t666 = 0.1e1 / (t855 * t892 + (pkin(3) * t845 + t725 * t764) * t781 + pkin(2) * t845);
t879 = t666 * t759;
t667 = 0.1e1 / (t853 * t891 + (pkin(3) * t843 + t726 * t764) * t783 + pkin(2) * t843);
t878 = t667 * t759;
t668 = 0.1e1 / (t851 * t890 + (pkin(3) * t841 + t727 * t764) * t785 + pkin(2) * t841);
t877 = t668 * t759;
t804 = 0.1e1 / pkin(3);
t876 = t673 * t804;
t875 = t674 * t804;
t874 = t675 * t804;
t873 = t676 * t804;
t734 = m(2) * rSges(2,2) - t889;
t835 = m(2) * rSges(2,1) + t790;
t872 = ((t812 * m(3) + t835) * t774 - t734 * t772) * t764;
t871 = ((t811 * m(3) + t835) * t782 - t734 * t776) * t764;
t870 = ((t810 * m(3) + t835) * t784 - t734 * t778) * t764;
t869 = ((t809 * m(3) + t835) * t786 - t734 * t780) * t764;
t802 = rSges(3,2) ^ 2;
t803 = rSges(3,1) ^ 2;
t733 = (t802 + t803) * m(3) + Icges(3,3);
t860 = t733 * t804;
t846 = t766 * t772;
t844 = t766 * t776;
t842 = t766 * t778;
t840 = t766 * t780;
t834 = pkin(2) * t859;
t833 = pkin(2) * t856;
t832 = pkin(2) * t854;
t831 = pkin(2) * t852;
t830 = -0.2e1 * pkin(2) * t898;
t829 = t665 * t872;
t828 = t666 * t871;
t827 = t667 * t870;
t826 = t668 * t869;
t735 = -rSges(3,2) * t889 + Icges(3,6);
t736 = rSges(3,1) * t889 - Icges(3,5);
t685 = t735 * t773 - t736 * t771;
t825 = t685 * t876;
t824 = t673 * t860;
t686 = t735 * t781 - t736 * t775;
t823 = t686 * t875;
t822 = t674 * t860;
t687 = t735 * t783 - t736 * t777;
t821 = t687 * t874;
t820 = t675 * t860;
t688 = t735 * t785 - t736 * t779;
t819 = t688 * t873;
t818 = t676 * t860;
t669 = t763 * t724 - t808 * t765;
t697 = t763 * t846 - t765 * t774;
t698 = t763 * t774 + t765 * t846;
t621 = -(t697 * t751 + t747 * t698) * t893 + (-t747 * t669 + t900 * t751) * t773 + t693 * t834;
t622 = (-t747 * t697 + t698 * t751) * t893 + (t669 * t751 + t900 * t747) * t773 - t689 * t834;
t791 = xP(4);
t756 = sin(t791);
t757 = cos(t791);
t794 = koppelP(4,2);
t798 = koppelP(4,1);
t715 = -t756 * t798 - t757 * t794;
t719 = -t756 * t794 + t757 * t798;
t590 = (t621 * t715 + t622 * t719) * t665;
t645 = -t689 * t857 - (t689 * t846 + t774 * t693) * t771;
t646 = -t693 * t857 - (-t774 * t689 + t693 * t846) * t771;
t605 = (t645 * t715 + t646 * t719) * t665;
t609 = (t653 * t715 + t654 * t719) * t876;
t557 = t590 * t759 + t605 * t872 + t609 * t897;
t670 = t763 * t728 - t807 * t765;
t699 = t763 * t844 - t765 * t782;
t702 = t763 * t782 + t765 * t844;
t623 = -(t699 * t752 + t748 * t702) * t892 + (-t748 * t670 + t901 * t752) * t781 + t694 * t833;
t626 = (-t748 * t699 + t702 * t752) * t892 + (t670 * t752 + t901 * t748) * t781 - t690 * t833;
t795 = koppelP(3,2);
t799 = koppelP(3,1);
t716 = -t756 * t799 - t757 * t795;
t720 = -t756 * t795 + t757 * t799;
t594 = (t623 * t716 + t626 * t720) * t666;
t647 = -t690 * t850 - (t690 * t844 + t782 * t694) * t775;
t650 = -t694 * t850 - (-t782 * t690 + t694 * t844) * t775;
t606 = (t647 * t716 + t650 * t720) * t666;
t610 = (t655 * t716 + t658 * t720) * t875;
t559 = t594 * t759 + t606 * t871 + t610 * t896;
t671 = t763 * t729 - t806 * t765;
t700 = t763 * t842 - t765 * t784;
t703 = t763 * t784 + t765 * t842;
t624 = -(t700 * t753 + t749 * t703) * t891 + (-t749 * t671 + t902 * t753) * t783 + t695 * t832;
t627 = (-t749 * t700 + t703 * t753) * t891 + (t671 * t753 + t902 * t749) * t783 - t691 * t832;
t796 = koppelP(2,2);
t800 = koppelP(2,1);
t717 = -t756 * t800 - t757 * t796;
t721 = -t756 * t796 + t757 * t800;
t595 = (t624 * t717 + t627 * t721) * t667;
t648 = -t691 * t849 - (t691 * t842 + t784 * t695) * t777;
t651 = -t695 * t849 - (-t784 * t691 + t695 * t842) * t777;
t607 = (t648 * t717 + t651 * t721) * t667;
t611 = (t656 * t717 + t659 * t721) * t874;
t560 = t595 * t759 + t607 * t870 + t611 * t895;
t672 = t763 * t730 - t805 * t765;
t701 = t763 * t840 - t765 * t786;
t704 = t763 * t786 + t765 * t840;
t625 = -(t701 * t754 + t750 * t704) * t890 + (-t750 * t672 + t903 * t754) * t785 + t696 * t831;
t628 = (-t750 * t701 + t704 * t754) * t890 + (t672 * t754 + t903 * t750) * t785 - t692 * t831;
t797 = koppelP(1,2);
t801 = koppelP(1,1);
t718 = -t756 * t801 - t757 * t797;
t722 = -t756 * t797 + t757 * t801;
t596 = (t625 * t718 + t628 * t722) * t668;
t649 = -t692 * t848 - (t692 * t840 + t786 * t696) * t779;
t652 = -t696 * t848 - (-t786 * t692 + t696 * t840) * t779;
t608 = (t649 * t718 + t652 * t722) * t668;
t612 = (t657 * t718 + t660 * t722) * t873;
t561 = t596 * t759 + t608 * t869 + t612 * t894;
t817 = t876 * t897;
t573 = t621 * t880 + t645 * t829 + t653 * t817;
t574 = t622 * t880 + t646 * t829 + t654 * t817;
t816 = t875 * t896;
t577 = t623 * t879 + t647 * t828 + t655 * t816;
t815 = t874 * t895;
t578 = t624 * t878 + t648 * t827 + t656 * t815;
t814 = t873 * t894;
t579 = t625 * t877 + t649 * t826 + t657 * t814;
t580 = t626 * t879 + t650 * t828 + t658 * t816;
t581 = t627 * t878 + t651 * t827 + t659 * t815;
t582 = t628 * t877 + t652 * t826 + t660 * t814;
t813 = Icges(3,1) + Icges(2,3) + (pkin(2) ^ 2 + t787 ^ 2 + t802) * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t793 = rSges(4,1);
t792 = rSges(4,2);
t755 = 0.2e1 * rSges(3,1) * t790;
t731 = (-t802 + t803) * m(3) + Icges(3,2) - Icges(3,1);
t708 = m(4) * (-t756 * t792 + t757 * t793);
t707 = m(4) * (-t756 * t793 - t757 * t792);
t664 = t731 * t762 + (t779 * t899 + t755) * t785 + t779 * t830 + t813;
t663 = t731 * t761 + (t777 * t899 + t755) * t783 + t777 * t830 + t813;
t662 = t731 * t760 + (t775 * t899 + t755) * t781 + t775 * t830 + t813;
t661 = t731 * t758 + (t771 * t899 + t755) * t773 + t771 * t830 + t813;
t588 = t660 * t818 + (t628 * t894 + t652 * t688) * t668;
t587 = t659 * t820 + (t627 * t895 + t651 * t687) * t667;
t586 = t658 * t822 + (t626 * t896 + t650 * t686) * t666;
t585 = t657 * t818 + (t625 * t894 + t649 * t688) * t668;
t584 = t656 * t820 + (t624 * t895 + t648 * t687) * t667;
t583 = t655 * t822 + (t623 * t896 + t647 * t686) * t666;
t576 = t654 * t824 + (t622 * t897 + t646 * t685) * t665;
t575 = t653 * t824 + (t621 * t897 + t645 * t685) * t665;
t572 = t660 * t819 + (t628 * t869 + t652 * t664) * t668;
t571 = t659 * t821 + (t627 * t870 + t651 * t663) * t667;
t570 = t658 * t823 + (t626 * t871 + t650 * t662) * t666;
t569 = t657 * t819 + (t625 * t869 + t649 * t664) * t668;
t568 = t656 * t821 + (t624 * t870 + t648 * t663) * t667;
t567 = t655 * t823 + (t623 * t871 + t647 * t662) * t666;
t566 = t654 * t825 + (t622 * t872 + t646 * t661) * t665;
t565 = t653 * t825 + (t621 * t872 + t645 * t661) * t665;
t564 = t596 * t894 + t608 * t688 + t612 * t733;
t563 = t595 * t895 + t607 * t687 + t611 * t733;
t562 = t594 * t896 + t606 * t686 + t610 * t733;
t558 = t590 * t897 + t605 * t685 + t609 * t733;
t556 = t596 * t869 + t608 * t664 + t612 * t688;
t555 = t595 * t870 + t607 * t663 + t611 * t687;
t554 = t594 * t871 + t606 * t662 + t610 * t686;
t553 = t590 * t872 + t605 * t661 + t609 * t685;
t552 = t582 + t581 + t580 + t574;
t551 = t579 + t578 + t577 + t573;
t550 = t561 + t560 + t559 + t557;
t1 = [m(4) + (t569 * t649 + t579 * t625) * t668 + (t568 * t648 + t578 * t624) * t667 + (t567 * t647 + t577 * t623) * t666 + (t565 * t645 + t573 * t621) * t665 + (t575 * t888 + t583 * t886 + t584 * t885 + t585 * t884) * t804, (t569 * t652 + t579 * t628) * t668 + (t568 * t651 + t578 * t627) * t667 + (t567 * t650 + t577 * t626) * t666 + (t565 * t646 + t573 * t622) * t665 + (t575 * t887 + t583 * t883 + t584 * t882 + t585 * t881) * t804, t551, t565 * t605 + t567 * t606 + t568 * t607 + t569 * t608 + t573 * t590 + t575 * t609 + t577 * t594 + t578 * t595 + t579 * t596 + t583 * t610 + t584 * t611 + t585 * t612 + t707; (t572 * t649 + t582 * t625) * t668 + (t571 * t648 + t581 * t624) * t667 + (t570 * t647 + t580 * t623) * t666 + (t566 * t645 + t574 * t621) * t665 + (t576 * t888 + t586 * t886 + t587 * t885 + t588 * t884) * t804, m(4) + (t572 * t652 + t582 * t628) * t668 + (t571 * t651 + t581 * t627) * t667 + (t570 * t650 + t580 * t626) * t666 + (t566 * t646 + t574 * t622) * t665 + (t576 * t887 + t586 * t883 + t587 * t882 + t588 * t881) * t804, t552, t566 * t605 + t570 * t606 + t571 * t607 + t572 * t608 + t574 * t590 + t576 * t609 + t580 * t594 + t581 * t595 + t582 * t596 + t586 * t610 + t587 * t611 + t588 * t612 + t708; t551, t552, 0.4e1 * m(1) + 0.4e1 * m(2) + 0.4e1 * m(3) + m(4), t550; t707 + (t556 * t649 + t561 * t625) * t668 + (t555 * t648 + t560 * t624) * t667 + (t554 * t647 + t559 * t623) * t666 + (t553 * t645 + t557 * t621) * t665 + (t558 * t888 + t562 * t886 + t563 * t885 + t564 * t884) * t804, t708 + (t556 * t652 + t561 * t628) * t668 + (t555 * t651 + t560 * t627) * t667 + (t554 * t650 + t559 * t626) * t666 + (t553 * t646 + t557 * t622) * t665 + (t558 * t887 + t562 * t883 + t563 * t882 + t564 * t881) * t804, t550, t561 * t596 + t556 * t608 + t564 * t612 + t560 * t595 + t555 * t607 + t563 * t611 + t559 * t594 + t554 * t606 + t562 * t610 + t557 * t590 + t553 * t605 + t558 * t609 + Icges(4,3) + m(4) * (t792 ^ 2 + t793 ^ 2);];
MX  = t1;
