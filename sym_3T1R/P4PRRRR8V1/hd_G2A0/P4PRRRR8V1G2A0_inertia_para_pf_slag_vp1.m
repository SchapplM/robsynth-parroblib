% Calculate inertia matrix for parallel robot
% P4PRRRR8V1G2A0
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
% MX [4x4]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:06
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRRR8V1G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:02:25
% EndTime: 2020-08-07 11:02:29
% DurationCPUTime: 3.64s
% Computational Cost: add. (9086->404), mult. (24847->729), div. (1920->9), fcn. (22856->38), ass. (0->323)
t928 = m(3) / 0.2e1;
t927 = Icges(3,2) / 0.2e1;
t752 = sin(qJ(2,4));
t754 = cos(qJ(2,4));
t753 = cos(qJ(3,4));
t914 = pkin(2) * t753;
t712 = -t754 * pkin(5) + t752 * t914;
t748 = sin(pkin(3));
t750 = cos(pkin(3));
t751 = sin(qJ(3,4));
t873 = t751 * t750;
t691 = pkin(2) * t873 + t712 * t748;
t926 = 0.1e1 / t691;
t760 = sin(qJ(2,3));
t766 = cos(qJ(2,3));
t765 = cos(qJ(3,3));
t913 = pkin(2) * t765;
t714 = -t766 * pkin(5) + t760 * t913;
t759 = sin(qJ(3,3));
t871 = t759 * t750;
t695 = pkin(2) * t871 + t714 * t748;
t925 = 0.1e1 / t695;
t762 = sin(qJ(2,2));
t768 = cos(qJ(2,2));
t767 = cos(qJ(3,2));
t912 = pkin(2) * t767;
t715 = -t768 * pkin(5) + t762 * t912;
t761 = sin(qJ(3,2));
t869 = t761 * t750;
t696 = pkin(2) * t869 + t715 * t748;
t924 = 0.1e1 / t696;
t764 = sin(qJ(2,1));
t770 = cos(qJ(2,1));
t769 = cos(qJ(3,1));
t911 = pkin(2) * t769;
t716 = -t770 * pkin(5) + t764 * t911;
t763 = sin(qJ(3,1));
t867 = t763 * t750;
t697 = pkin(2) * t867 + t716 * t748;
t923 = 0.1e1 / t697;
t922 = m(3) * rSges(3,3);
t789 = rSges(3,2) ^ 2;
t790 = rSges(3,1) ^ 2;
t921 = (-t789 + t790) * t928 - Icges(3,1) / 0.2e1 + t927;
t820 = rSges(3,1) * t753 - rSges(3,2) * t751;
t678 = t820 * t750 - t748 * t752 * (t751 * rSges(3,1) + t753 * rSges(3,2));
t920 = m(3) * t678;
t819 = rSges(3,1) * t765 - rSges(3,2) * t759;
t679 = t819 * t750 - t748 * t760 * (t759 * rSges(3,1) + t765 * rSges(3,2));
t919 = m(3) * t679;
t818 = rSges(3,1) * t767 - rSges(3,2) * t761;
t680 = t818 * t750 - t748 * t762 * (t761 * rSges(3,1) + t767 * rSges(3,2));
t918 = m(3) * t680;
t817 = rSges(3,1) * t769 - rSges(3,2) * t763;
t681 = t817 * t750 - t748 * t764 * (t763 * rSges(3,1) + t769 * rSges(3,2));
t917 = m(3) * t681;
t791 = 0.1e1 / pkin(2);
t916 = m(3) * t791;
t915 = pkin(2) * t748;
t713 = pkin(5) * t752 + t754 * t914;
t747 = sin(pkin(6));
t749 = cos(pkin(6));
t815 = -t712 * t750 + t751 * t915;
t654 = -t747 * t713 + t815 * t749;
t755 = legFrame(4,2);
t732 = sin(t755);
t736 = cos(t755);
t642 = t654 * t732 + t736 * t691;
t910 = t642 * t926;
t643 = -t654 * t736 + t732 * t691;
t909 = t643 * t926;
t717 = pkin(5) * t760 + t766 * t913;
t814 = -t714 * t750 + t759 * t915;
t656 = -t747 * t717 + t814 * t749;
t756 = legFrame(3,2);
t733 = sin(t756);
t737 = cos(t756);
t644 = t656 * t733 + t737 * t695;
t908 = t644 * t925;
t645 = -t656 * t737 + t733 * t695;
t907 = t645 * t925;
t718 = pkin(5) * t762 + t768 * t912;
t813 = -t715 * t750 + t761 * t915;
t657 = -t747 * t718 + t813 * t749;
t757 = legFrame(2,2);
t734 = sin(t757);
t738 = cos(t757);
t646 = t657 * t734 + t738 * t696;
t906 = t646 * t924;
t647 = -t657 * t738 + t734 * t696;
t905 = t647 * t924;
t719 = pkin(5) * t764 + t770 * t911;
t812 = -t716 * t750 + t763 * t915;
t658 = -t747 * t719 + t812 * t749;
t758 = legFrame(1,2);
t735 = sin(t758);
t739 = cos(t758);
t648 = t658 * t735 + t739 * t697;
t904 = t648 * t923;
t649 = -t658 * t739 + t735 * t697;
t903 = t649 * t923;
t655 = t749 * t713 + t815 * t747;
t902 = t655 * t926;
t659 = t749 * t717 + t814 * t747;
t901 = t659 * t925;
t660 = t749 * t718 + t813 * t747;
t900 = t660 * t924;
t661 = t749 * t719 + t812 * t747;
t899 = t661 * t923;
t728 = m(2) * rSges(2,2) - t922;
t773 = m(2) * rSges(2,1);
t898 = ((t820 * m(3) + t773) * t754 - t728 * t752) * t748;
t897 = ((t819 * m(3) + t773) * t766 - t728 * t760) * t748;
t896 = ((t818 * m(3) + t773) * t768 - t728 * t762) * t748;
t895 = ((t817 * m(3) + t773) * t770 - t728 * t764) * t748;
t894 = t926 / t753;
t893 = t925 / t765;
t892 = t924 / t767;
t891 = t923 / t769;
t743 = m(1) + m(2) + m(3);
t890 = t926 * t743;
t889 = t925 * t743;
t888 = t924 * t743;
t887 = t923 * t743;
t729 = -rSges(3,2) * t922 + Icges(3,6);
t730 = rSges(3,1) * t922 - Icges(3,5);
t698 = t729 * t753 - t730 * t751;
t886 = t698 * t791;
t699 = t729 * t765 - t730 * t759;
t885 = t699 * t791;
t700 = t729 * t767 - t730 * t761;
t884 = t700 * t791;
t701 = t729 * t769 - t730 * t763;
t883 = t701 * t791;
t865 = t789 + t790;
t723 = t865 * m(3) + Icges(3,3);
t882 = t723 * t791;
t881 = t750 * t752;
t880 = t750 * t754;
t879 = t750 * t760;
t878 = t750 * t762;
t877 = t750 * t764;
t876 = t750 * t766;
t875 = t750 * t768;
t874 = t750 * t770;
t872 = t751 * t754;
t870 = t759 * t766;
t868 = t761 * t768;
t866 = t763 * t770;
t864 = t926 * t920;
t863 = t678 * t916;
t862 = t925 * t919;
t861 = t679 * t916;
t860 = t924 * t918;
t859 = t680 * t916;
t858 = t923 * t917;
t857 = t681 * t916;
t662 = -(-t747 * t752 + t749 * t880) * t914 - pkin(5) * (t747 * t754 + t749 * t881);
t856 = t662 * t894;
t664 = -(-t747 * t760 + t749 * t876) * t913 - pkin(5) * (t747 * t766 + t749 * t879);
t855 = t664 * t893;
t665 = -(-t747 * t762 + t749 * t875) * t912 - pkin(5) * (t747 * t768 + t749 * t878);
t854 = t665 * t892;
t666 = -(-t747 * t764 + t749 * t874) * t911 - pkin(5) * (t747 * t770 + t749 * t877);
t853 = t666 * t891;
t803 = t748 * t753 + t752 * t873;
t671 = -t747 * t872 - t803 * t749;
t852 = t671 * t894;
t802 = t748 * t765 + t760 * t871;
t675 = -t747 * t870 - t802 * t749;
t851 = t675 * t893;
t801 = t748 * t767 + t762 * t869;
t676 = -t747 * t868 - t801 * t749;
t850 = t676 * t892;
t800 = t748 * t769 + t764 * t867;
t677 = -t747 * t866 - t800 * t749;
t849 = t677 * t891;
t848 = t926 * t898;
t847 = t925 * t897;
t846 = t924 * t896;
t845 = t923 * t895;
t844 = t732 * t894;
t843 = t736 * t894;
t842 = t733 * t893;
t841 = t737 * t893;
t840 = t734 * t892;
t839 = t738 * t892;
t838 = t735 * t891;
t837 = t739 * t891;
t663 = (t747 * t880 + t749 * t752) * t914 + (t747 * t881 - t749 * t754) * pkin(5);
t836 = t663 * t844;
t835 = t663 * t843;
t667 = (t747 * t876 + t749 * t760) * t913 + (t747 * t879 - t749 * t766) * pkin(5);
t834 = t667 * t842;
t833 = t667 * t841;
t668 = (t747 * t875 + t749 * t762) * t912 + (t747 * t878 - t749 * t768) * pkin(5);
t832 = t668 * t840;
t831 = t668 * t839;
t669 = (t747 * t874 + t749 * t764) * t911 + (t747 * t877 - t749 * t770) * pkin(5);
t830 = t669 * t838;
t829 = t669 * t837;
t670 = t803 * t747 - t749 * t872;
t828 = t670 * t844;
t827 = t670 * t843;
t672 = t802 * t747 - t749 * t870;
t826 = t672 * t842;
t825 = t672 * t841;
t673 = t801 * t747 - t749 * t868;
t824 = t673 * t840;
t823 = t673 * t839;
t674 = t800 * t747 - t749 * t866;
t822 = t674 * t838;
t821 = t674 * t837;
t816 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + (0.2e1 * rSges(3,3) ^ 2 + t865) * t928 + t927 + Icges(3,1) / 0.2e1;
t731 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t777 = 0.2e1 * qJ(3,4);
t650 = cos(t777) * t921 + t731 * sin(t777) + t816;
t811 = t650 * t670 + t663 * t886;
t778 = 0.2e1 * qJ(3,3);
t651 = cos(t778) * t921 + t731 * sin(t778) + t816;
t810 = t651 * t672 + t667 * t885;
t779 = 0.2e1 * qJ(3,2);
t652 = cos(t779) * t921 + t731 * sin(t779) + t816;
t809 = t652 * t673 + t668 * t884;
t780 = 0.2e1 * qJ(3,1);
t653 = cos(t780) * t921 + t731 * sin(t780) + t816;
t808 = t653 * t674 + t669 * t883;
t807 = t663 * t882 + t670 * t698;
t806 = t667 * t882 + t672 * t699;
t805 = t668 * t882 + t673 * t700;
t804 = t669 * t882 + t674 * t701;
t774 = xP(4);
t740 = sin(t774);
t741 = cos(t774);
t781 = koppelP(4,2);
t785 = koppelP(4,1);
t704 = -t740 * t785 - t741 * t781;
t708 = -t740 * t781 + t741 * t785;
t799 = (-t704 * t736 + t708 * t732) * t894;
t782 = koppelP(3,2);
t786 = koppelP(3,1);
t705 = -t740 * t786 - t741 * t782;
t709 = -t740 * t782 + t741 * t786;
t798 = (-t705 * t737 + t709 * t733) * t893;
t783 = koppelP(2,2);
t787 = koppelP(2,1);
t706 = -t740 * t787 - t741 * t783;
t710 = -t740 * t783 + t741 * t787;
t797 = (-t706 * t738 + t710 * t734) * t892;
t784 = koppelP(1,2);
t788 = koppelP(1,1);
t707 = -t740 * t788 - t741 * t784;
t711 = -t740 * t784 + t741 * t788;
t796 = (-t707 * t739 + t711 * t735) * t891;
t795 = t663 * t863 + t670 * t898;
t794 = t667 * t861 + t672 * t897;
t793 = t668 * t859 + t673 * t896;
t792 = t669 * t857 + t674 * t895;
t776 = rSges(4,1);
t775 = rSges(4,2);
t703 = m(4) * (-t740 * t775 + t741 * t776);
t702 = m(4) * (-t740 * t776 - t741 * t775);
t641 = t674 * t796;
t640 = t673 * t797;
t639 = t672 * t798;
t638 = t670 * t799;
t637 = t791 * t669 * t796;
t636 = t791 * t668 * t797;
t635 = t791 * t667 * t798;
t634 = t791 * t663 * t799;
t633 = (t648 * t711 + t649 * t707) * t923;
t632 = (t646 * t710 + t647 * t706) * t924;
t631 = (t644 * t709 + t645 * t705) * t925;
t630 = (t642 * t708 + t643 * t704) * t926;
t629 = t661 * t858 + (t666 * t882 + t677 * t701) * t891;
t628 = t660 * t860 + (t665 * t882 + t676 * t700) * t892;
t627 = t659 * t862 + (t664 * t882 + t675 * t699) * t893;
t626 = t655 * t864 + (t662 * t882 + t671 * t698) * t894;
t625 = t661 * t887 + (t666 * t857 + t677 * t895) * t891;
t624 = t660 * t888 + (t665 * t859 + t676 * t896) * t892;
t623 = t659 * t889 + (t664 * t861 + t675 * t897) * t893;
t622 = t655 * t890 + (t662 * t863 + t671 * t898) * t894;
t621 = t661 * t845 + (t653 * t677 + t666 * t883) * t891;
t620 = t660 * t846 + (t652 * t676 + t665 * t884) * t892;
t619 = t659 * t847 + (t651 * t675 + t664 * t885) * t893;
t618 = t649 * t858 - t804 * t837;
t617 = t648 * t858 + t804 * t838;
t616 = t647 * t860 - t805 * t839;
t615 = t646 * t860 + t805 * t840;
t614 = t645 * t862 - t806 * t841;
t613 = t644 * t862 + t806 * t842;
t612 = t655 * t848 + (t650 * t671 + t662 * t886) * t894;
t611 = t649 * t887 - t792 * t837;
t610 = t648 * t887 + t792 * t838;
t609 = t647 * t888 - t793 * t839;
t608 = t646 * t888 + t793 * t840;
t607 = t645 * t889 - t794 * t841;
t606 = t644 * t889 + t794 * t842;
t605 = t643 * t864 - t807 * t843;
t604 = t642 * t864 + t807 * t844;
t603 = t643 * t890 - t795 * t843;
t602 = t642 * t890 + t795 * t844;
t601 = t649 * t845 - t808 * t837;
t600 = t648 * t845 + t808 * t838;
t599 = t647 * t846 - t809 * t839;
t598 = t646 * t846 + t809 * t840;
t597 = t645 * t847 - t810 * t841;
t596 = t644 * t847 + t810 * t842;
t595 = t643 * t848 - t811 * t843;
t594 = t642 * t848 + t811 * t844;
t593 = t633 * t917 + t637 * t723 + t641 * t701;
t592 = t632 * t918 + t636 * t723 + t640 * t700;
t591 = t631 * t919 + t635 * t723 + t639 * t699;
t590 = t633 * t743 + t637 * t917 + t641 * t895;
t589 = t632 * t743 + t636 * t918 + t640 * t896;
t588 = t631 * t743 + t635 * t919 + t639 * t897;
t587 = t630 * t920 + t634 * t723 + t638 * t698;
t586 = t630 * t743 + t634 * t920 + t638 * t898;
t585 = t633 * t895 + t637 * t701 + t641 * t653;
t584 = t632 * t896 + t636 * t700 + t640 * t652;
t583 = t631 * t897 + t635 * t699 + t639 * t651;
t582 = t630 * t898 + t634 * t698 + t638 * t650;
t1 = [t611 * t903 - t601 * t821 + t609 * t905 - t599 * t823 + t607 * t907 - t597 * t825 + t603 * t909 - t595 * t827 + m(4) + (-t605 * t835 - t614 * t833 - t616 * t831 - t618 * t829) * t791, t611 * t904 + t601 * t822 + t609 * t906 + t599 * t824 + t607 * t908 + t597 * t826 + t603 * t910 + t595 * t828 + (t605 * t836 + t614 * t834 + t616 * t832 + t618 * t830) * t791, t595 * t852 + t597 * t851 + t599 * t850 + t601 * t849 + t603 * t902 + t607 * t901 + t609 * t900 + t611 * t899 + (t605 * t856 + t614 * t855 + t616 * t854 + t618 * t853) * t791, t595 * t638 + t597 * t639 + t599 * t640 + t601 * t641 + t603 * t630 + t605 * t634 + t607 * t631 + t609 * t632 + t611 * t633 + t614 * t635 + t616 * t636 + t618 * t637 + t702; t610 * t903 - t600 * t821 + t608 * t905 - t598 * t823 + t606 * t907 - t596 * t825 + t602 * t909 - t594 * t827 + (-t604 * t835 - t613 * t833 - t615 * t831 - t617 * t829) * t791, t610 * t904 + t600 * t822 + t608 * t906 + t598 * t824 + t606 * t908 + t596 * t826 + t602 * t910 + t594 * t828 + m(4) + (t604 * t836 + t613 * t834 + t615 * t832 + t617 * t830) * t791, t594 * t852 + t596 * t851 + t598 * t850 + t600 * t849 + t602 * t902 + t606 * t901 + t608 * t900 + t610 * t899 + (t604 * t856 + t613 * t855 + t615 * t854 + t617 * t853) * t791, t594 * t638 + t596 * t639 + t598 * t640 + t600 * t641 + t602 * t630 + t604 * t634 + t606 * t631 + t608 * t632 + t610 * t633 + t613 * t635 + t615 * t636 + t617 * t637 + t703; t625 * t903 - t621 * t821 + t624 * t905 - t620 * t823 + t623 * t907 - t619 * t825 + t622 * t909 - t612 * t827 + (-t626 * t835 - t627 * t833 - t628 * t831 - t629 * t829) * t791, t625 * t904 + t621 * t822 + t624 * t906 + t620 * t824 + t623 * t908 + t619 * t826 + t622 * t910 + t612 * t828 + (t626 * t836 + t627 * t834 + t628 * t832 + t629 * t830) * t791, t612 * t852 + t619 * t851 + t620 * t850 + t621 * t849 + t622 * t902 + t623 * t901 + t624 * t900 + t625 * t899 + m(4) + (t626 * t856 + t627 * t855 + t628 * t854 + t629 * t853) * t791, t612 * t638 + t619 * t639 + t620 * t640 + t621 * t641 + t622 * t630 + t623 * t631 + t624 * t632 + t625 * t633 + t626 * t634 + t627 * t635 + t628 * t636 + t629 * t637; t590 * t903 - t585 * t821 + t589 * t905 - t584 * t823 + t588 * t907 - t583 * t825 + t586 * t909 - t582 * t827 + t702 + (-t587 * t835 - t591 * t833 - t592 * t831 - t593 * t829) * t791, t590 * t904 + t585 * t822 + t589 * t906 + t584 * t824 + t588 * t908 + t583 * t826 + t586 * t910 + t582 * t828 + t703 + (t587 * t836 + t591 * t834 + t592 * t832 + t593 * t830) * t791, t582 * t852 + t583 * t851 + t584 * t850 + t585 * t849 + t586 * t902 + t588 * t901 + t589 * t900 + t590 * t899 + (t587 * t856 + t591 * t855 + t592 * t854 + t593 * t853) * t791, t590 * t633 + t585 * t641 + t593 * t637 + t589 * t632 + t584 * t640 + t592 * t636 + t588 * t631 + t583 * t639 + t591 * t635 + t586 * t630 + t582 * t638 + t587 * t634 + Icges(4,3) + m(4) * (t775 ^ 2 + t776 ^ 2);];
MX  = t1;
