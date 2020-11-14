% Calculate vector of centrifugal and coriolis load on the joints for
% P4PRRRR1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% taucX [4x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 19:06
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:00:51
% EndTime: 2020-03-02 19:00:59
% DurationCPUTime: 8.66s
% Computational Cost: add. (3948->336), mult. (8530->725), div. (3852->18), fcn. (6884->34), ass. (0->291)
t744 = sin(qJ(3,1));
t750 = cos(qJ(3,1));
t903 = (rSges(3,1) * t750 - rSges(3,2) * t744) * m(3);
t742 = sin(qJ(3,2));
t748 = cos(qJ(3,2));
t902 = (rSges(3,1) * t748 - rSges(3,2) * t742) * m(3);
t740 = sin(qJ(3,3));
t746 = cos(qJ(3,3));
t901 = (rSges(3,1) * t746 - rSges(3,2) * t740) * m(3);
t732 = sin(qJ(3,4));
t734 = cos(qJ(3,4));
t900 = (rSges(3,1) * t734 - rSges(3,2) * t732) * m(3);
t710 = t734 ^ 2;
t711 = 0.1e1 / t734;
t713 = t711 / t710;
t899 = t713 * t732;
t718 = t746 ^ 2;
t719 = 0.1e1 / t746;
t721 = t719 / t718;
t898 = t721 * t740;
t722 = t748 ^ 2;
t723 = 0.1e1 / t748;
t725 = t723 / t722;
t897 = t725 * t742;
t726 = t750 ^ 2;
t727 = 0.1e1 / t750;
t729 = t727 / t726;
t896 = t729 * t744;
t759 = xP(4);
t707 = sin(t759);
t708 = cos(t759);
t766 = koppelP(4,2);
t770 = koppelP(4,1);
t670 = t707 * t770 + t708 * t766;
t674 = -t707 * t766 + t708 * t770;
t736 = legFrame(4,2);
t699 = sin(t736);
t703 = cos(t736);
t733 = sin(qJ(2,4));
t709 = 0.1e1 / t733;
t712 = 0.1e1 / t734 ^ 2;
t735 = cos(qJ(2,4));
t755 = xDP(4);
t757 = xDP(2);
t758 = xDP(1);
t776 = 0.1e1 / pkin(2);
t756 = xDP(3);
t846 = t732 * t756;
t623 = (-t735 * t712 * t846 + (-t703 * (-t670 * t755 + t758) + t699 * (t674 * t755 + t757)) * t711) * t776 * t709;
t622 = t623 ^ 2;
t839 = t756 * t776;
t804 = 0.2e1 * m(3) * t839;
t694 = rSges(3,2) * t804;
t886 = rSges(3,3) * m(3);
t695 = m(2) * rSges(2,2) - t886;
t754 = m(2) * rSges(2,1);
t791 = rSges(3,1) * t804;
t731 = t756 ^ 2;
t777 = 0.1e1 / pkin(2) ^ 2;
t848 = t731 * t777;
t873 = t711 * t732;
t881 = t623 * t735;
t578 = (-t622 * t754 - (t712 * t848 + t622) * t900) * t733 - (t623 * t695 + t791 * t873 + t694) * t881;
t822 = t711 * t839;
t847 = t732 * t734;
t871 = t712 * t776;
t590 = ((t710 * pkin(2) * t881 - t733 * t846) * t623 * t871 + (-t733 * t623 * t847 + t735 * t822) * t713 * t839) * t709;
t849 = t731 * t776;
t614 = (-pkin(2) * t622 * t734 - t713 * t849) * t709;
t654 = (t754 + t900) * t735 - t733 * t695;
t714 = m(1) + m(2) + m(3);
t795 = t848 * t899;
t836 = m(3) * (t732 * rSges(3,1) + t734 * rSges(3,2)) * t733;
t582 = -t654 * t590 - t714 * t614 - t795 * t836;
t895 = t582 + t578;
t767 = koppelP(3,2);
t771 = koppelP(3,1);
t671 = t707 * t771 + t708 * t767;
t675 = -t707 * t767 + t708 * t771;
t737 = legFrame(3,2);
t700 = sin(t737);
t704 = cos(t737);
t741 = sin(qJ(2,3));
t715 = 0.1e1 / t741;
t720 = 0.1e1 / t746 ^ 2;
t747 = cos(qJ(2,3));
t844 = t740 * t756;
t627 = (-t747 * t720 * t844 + (-t704 * (-t671 * t755 + t758) + t700 * (t675 * t755 + t757)) * t719) * t776 * t715;
t624 = t627 ^ 2;
t858 = t719 * t740;
t880 = t627 * t747;
t579 = (-t624 * t754 - (t720 * t848 + t624) * t901) * t741 - (t627 * t695 + t791 * t858 + t694) * t880;
t814 = t719 * t839;
t845 = t740 * t746;
t856 = t720 * t776;
t593 = ((t718 * pkin(2) * t880 - t741 * t844) * t627 * t856 + (-t741 * t627 * t845 + t747 * t814) * t721 * t839) * t715;
t615 = (-pkin(2) * t624 * t746 - t721 * t849) * t715;
t655 = (t754 + t901) * t747 - t741 * t695;
t794 = t848 * t898;
t835 = m(3) * (t740 * rSges(3,1) + t746 * rSges(3,2)) * t741;
t583 = -t655 * t593 - t714 * t615 - t794 * t835;
t894 = t583 + t579;
t768 = koppelP(2,2);
t772 = koppelP(2,1);
t672 = t707 * t772 + t708 * t768;
t676 = -t707 * t768 + t708 * t772;
t738 = legFrame(2,2);
t701 = sin(t738);
t705 = cos(t738);
t743 = sin(qJ(2,2));
t716 = 0.1e1 / t743;
t724 = 0.1e1 / t748 ^ 2;
t749 = cos(qJ(2,2));
t842 = t742 * t756;
t628 = (-t749 * t724 * t842 + (-t705 * (-t672 * t755 + t758) + t701 * (t676 * t755 + t757)) * t723) * t776 * t716;
t625 = t628 ^ 2;
t855 = t723 * t742;
t879 = t628 * t749;
t580 = (-t625 * t754 - (t724 * t848 + t625) * t902) * t743 - (t628 * t695 + t791 * t855 + t694) * t879;
t812 = t723 * t839;
t843 = t742 * t748;
t853 = t724 * t776;
t591 = ((t722 * pkin(2) * t879 - t743 * t842) * t628 * t853 + (-t743 * t628 * t843 + t749 * t812) * t725 * t839) * t716;
t616 = (-pkin(2) * t625 * t748 - t725 * t849) * t716;
t656 = (t754 + t902) * t749 - t743 * t695;
t793 = t848 * t897;
t834 = m(3) * (t742 * rSges(3,1) + t748 * rSges(3,2)) * t743;
t584 = -t656 * t591 - t714 * t616 - t793 * t834;
t893 = t584 + t580;
t769 = koppelP(1,2);
t773 = koppelP(1,1);
t673 = t707 * t773 + t708 * t769;
t677 = -t707 * t769 + t708 * t773;
t739 = legFrame(1,2);
t702 = sin(t739);
t706 = cos(t739);
t745 = sin(qJ(2,1));
t717 = 0.1e1 / t745;
t728 = 0.1e1 / t750 ^ 2;
t751 = cos(qJ(2,1));
t840 = t744 * t756;
t629 = (-t751 * t728 * t840 + (-t706 * (-t673 * t755 + t758) + t702 * (t677 * t755 + t757)) * t727) * t776 * t717;
t626 = t629 ^ 2;
t852 = t727 * t744;
t878 = t629 * t751;
t581 = (-t626 * t754 - (t728 * t848 + t626) * t903) * t745 - (t629 * t695 + t791 * t852 + t694) * t878;
t810 = t727 * t839;
t841 = t744 * t750;
t850 = t728 * t776;
t592 = ((t726 * pkin(2) * t878 - t745 * t840) * t629 * t850 + (-t745 * t629 * t841 + t751 * t810) * t729 * t839) * t717;
t617 = (-pkin(2) * t626 * t750 - t729 * t849) * t717;
t657 = (t754 + t903) * t751 - t745 * t695;
t792 = t848 * t896;
t833 = m(3) * (t744 * rSges(3,1) + t750 * rSges(3,2)) * t745;
t585 = -t657 * t592 - t714 * t617 - t792 * t833;
t892 = t585 + t581;
t698 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t765 = 0.2e1 * qJ(3,1);
t774 = rSges(3,2) ^ 2;
t775 = rSges(3,1) ^ 2;
t837 = t774 + t775;
t782 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + (0.2e1 * rSges(3,3) ^ 2 + t837) * m(3) / 0.2e1;
t688 = (-t774 + t775) * m(3) - Icges(3,1) + Icges(3,2);
t885 = t688 / 0.2e1;
t653 = cos(t765) * t885 - t698 * sin(t765) + t782;
t696 = rSges(3,2) * t886 - Icges(3,6);
t697 = rSges(3,1) * t886 - Icges(3,5);
t661 = -t696 * t750 - t697 * t744;
t577 = -t653 * t592 - t657 * t617 + t661 * t792;
t805 = t688 * t841;
t883 = t697 / 0.4e1;
t884 = -t696 / 0.4e1;
t589 = (t744 * t884 + t750 * t883) * t810 + (t805 / 0.2e1 + (t726 - 0.1e1 / 0.2e1) * t698) * t629;
t838 = t756 * t777;
t796 = t589 * t717 * t838;
t851 = t727 * t776;
t815 = t717 * t851;
t891 = t577 * t815 - 0.4e1 * t728 * t796;
t764 = 0.2e1 * qJ(3,2);
t652 = cos(t764) * t885 - t698 * sin(t764) + t782;
t660 = -t696 * t748 - t697 * t742;
t576 = -t652 * t591 - t656 * t616 + t660 * t793;
t806 = t688 * t843;
t588 = (t742 * t884 + t748 * t883) * t812 + (t806 / 0.2e1 + (t722 - 0.1e1 / 0.2e1) * t698) * t628;
t797 = t588 * t716 * t838;
t854 = t723 * t776;
t817 = t716 * t854;
t890 = t576 * t817 - 0.4e1 * t724 * t797;
t763 = 0.2e1 * qJ(3,3);
t651 = cos(t763) * t885 - t698 * sin(t763) + t782;
t659 = -t696 * t746 - t697 * t740;
t575 = -t651 * t593 - t655 * t615 + t659 * t794;
t807 = t688 * t845;
t587 = (t740 * t884 + t746 * t883) * t814 + (t807 / 0.2e1 + (t718 - 0.1e1 / 0.2e1) * t698) * t627;
t798 = t587 * t715 * t838;
t857 = t719 * t776;
t819 = t715 * t857;
t889 = t575 * t819 - 0.4e1 * t720 * t798;
t762 = 0.2e1 * qJ(3,4);
t650 = cos(t762) * t885 - t698 * sin(t762) + t782;
t658 = -t696 * t734 - t697 * t732;
t574 = -t650 * t590 - t654 * t614 + t658 * t795;
t808 = t688 * t847;
t586 = (t732 * t884 + t734 * t883) * t822 + (t808 / 0.2e1 + (t710 - 0.1e1 / 0.2e1) * t698) * t623;
t799 = t586 * t709 * t838;
t872 = t711 * t776;
t823 = t709 * t872;
t888 = t574 * t823 - 0.4e1 * t712 * t799;
t730 = t755 ^ 2;
t887 = 0.2e1 * t698;
t882 = m(4) * t730;
t877 = t709 * t578;
t876 = t709 * t582;
t875 = t709 * t732;
t874 = t709 * t730;
t870 = t715 * t579;
t869 = t715 * t583;
t868 = t715 * t740;
t867 = t715 * t730;
t866 = t716 * t580;
t865 = t716 * t584;
t864 = t716 * t742;
t863 = t716 * t730;
t862 = t717 * t581;
t861 = t717 * t585;
t860 = t717 * t744;
t859 = t717 * t730;
t832 = t699 * t872;
t831 = t700 * t857;
t830 = t701 * t854;
t829 = t702 * t851;
t828 = t703 * t872;
t827 = t704 * t857;
t826 = t705 * t854;
t825 = t706 * t851;
t821 = t735 * t871;
t813 = t747 * t856;
t811 = t749 * t853;
t809 = t751 * t850;
t761 = rSges(4,1);
t760 = rSges(4,2);
t693 = -t837 * m(3) - Icges(3,3);
t669 = t745 * t702 + t706 * t751;
t668 = -t702 * t751 + t706 * t745;
t667 = t743 * t701 + t705 * t749;
t666 = -t701 * t749 + t705 * t743;
t665 = t741 * t700 + t704 * t747;
t664 = -t700 * t747 + t704 * t741;
t663 = t733 * t699 + t703 * t735;
t662 = -t699 * t735 + t703 * t733;
t649 = (t706 * t673 + t702 * t677) * t815;
t648 = (t705 * t672 + t701 * t676) * t817;
t647 = (t704 * t671 + t700 * t675) * t819;
t646 = (t703 * t670 + t699 * t674) * t823;
t645 = (-t657 * t825 + t669 * t714) * t717;
t644 = (t657 * t829 + t668 * t714) * t717;
t643 = (-t656 * t826 + t667 * t714) * t716;
t642 = (t656 * t830 + t666 * t714) * t716;
t641 = (-t655 * t827 + t665 * t714) * t715;
t640 = (t655 * t831 + t664 * t714) * t715;
t639 = (-t654 * t828 + t663 * t714) * t709;
t638 = (t654 * t832 + t662 * t714) * t709;
t637 = (t668 * t677 - t669 * t673) * t717;
t636 = (t666 * t676 - t667 * t672) * t716;
t635 = (t664 * t675 - t665 * t671) * t715;
t634 = (t662 * t674 - t663 * t670) * t709;
t633 = -t833 * t851 + (-t657 * t809 + t714 * t727) * t860;
t632 = -t834 * t854 + (-t656 * t811 + t714 * t723) * t864;
t631 = -t835 * t857 + (-t655 * t813 + t714 * t719) * t868;
t630 = -t836 * t872 + (-t654 * t821 + t711 * t714) * t875;
t613 = (-t653 * t825 + t657 * t669) * t717;
t612 = (t653 * t829 + t657 * t668) * t717;
t611 = (-t652 * t826 + t656 * t667) * t716;
t610 = (t652 * t830 + t656 * t666) * t716;
t609 = (-t651 * t827 + t655 * t665) * t715;
t608 = (t651 * t831 + t655 * t664) * t715;
t607 = (-t650 * t828 + t654 * t663) * t709;
t606 = (t650 * t832 + t654 * t662) * t709;
t605 = t661 * t851 + (-t653 * t809 + t657 * t727) * t860;
t604 = t660 * t854 + (-t652 * t811 + t656 * t723) * t864;
t603 = t659 * t857 + (-t651 * t813 + t655 * t719) * t868;
t602 = t658 * t872 + (-t650 * t821 + t654 * t711) * t875;
t601 = t637 * t714 + t649 * t657;
t600 = t636 * t714 + t648 * t656;
t599 = t635 * t714 + t647 * t655;
t598 = t634 * t714 + t646 * t654;
t597 = t637 * t657 + t649 * t653;
t596 = t636 * t656 + t648 * t652;
t595 = t635 * t655 + t647 * t651;
t594 = t634 * t654 + t646 * t650;
t1 = [(-(-t613 * t825 + t645 * t669) * t677 - (t613 * t829 + t645 * t668) * t673) * t859 + t669 * t861 + t669 * t862 + (-(-t611 * t826 + t643 * t667) * t676 - (t611 * t830 + t643 * t666) * t672) * t863 + t667 * t865 + t667 * t866 + (-(-t609 * t827 + t641 * t665) * t675 - (t609 * t831 + t641 * t664) * t671) * t867 + t665 * t869 + t665 * t870 + (-(-t607 * t828 + t639 * t663) * t674 - (t607 * t832 + t639 * t662) * t670) * t874 + t663 * t876 + t663 * t877 + (t707 * t760 - t708 * t761) * t882 - t891 * t706 - t890 * t705 - t889 * t704 - t888 * t703; (-(-t612 * t825 + t644 * t669) * t677 - (t612 * t829 + t644 * t668) * t673) * t859 + t668 * t861 + t668 * t862 + (-(-t610 * t826 + t642 * t667) * t676 - (t610 * t830 + t642 * t666) * t672) * t863 + t666 * t865 + t666 * t866 + (-(-t608 * t827 + t640 * t665) * t675 - (t608 * t831 + t640 * t664) * t671) * t867 + t664 * t869 + t664 * t870 + (-(-t606 * t828 + t638 * t663) * t674 - (t606 * t832 + t638 * t662) * t670) * t874 + t662 * t876 + t662 * t877 - (t707 * t761 + t708 * t760) * t882 + t891 * t702 + t890 * t701 + t889 * t700 + t888 * t699; (-(-t602 * t828 + t630 * t663) * t674 - (t602 * t832 + t630 * t662) * t670) * t874 + (-(-t605 * t825 + t633 * t669) * t677 - (t605 * t829 + t633 * t668) * t673) * t859 + (-(-t604 * t826 + t632 * t667) * t676 - (t604 * t830 + t632 * t666) * t672) * t863 + (-(-t603 * t827 + t631 * t665) * t675 - (t603 * t831 + t631 * t664) * t671) * t867 + 0.4e1 * t735 * t799 * t899 + 0.4e1 * t751 * t796 * t896 + 0.4e1 * t749 * t797 * t897 + 0.4e1 * t747 * t798 * t898 - t577 * t809 * t860 - t576 * t811 * t864 - t575 * t813 * t868 - t574 * t821 * t875 + (-t658 * t590 + t614 * t836 - t693 * t795 + (t710 * t887 - t698 + t808) * t622) * t872 + (-t659 * t593 + t615 * t835 - t693 * t794 + (t718 * t887 - t698 + t807) * t624) * t857 + (-t660 * t591 + t616 * t834 - t693 * t793 + (t722 * t887 - t698 + t806) * t625) * t854 + (-t661 * t592 + t617 * t833 - t693 * t792 + (t726 * t887 - t698 + t805) * t626) * t851 + t895 * t709 * t873 + t894 * t715 * t858 + t893 * t716 * t855 + t892 * t717 * t852; (-(-t597 * t825 + t601 * t669) * t677 - (t597 * t829 + t601 * t668) * t673) * t859 + (-(-t596 * t826 + t600 * t667) * t676 - (t596 * t830 + t600 * t666) * t672) * t863 + (-(-t595 * t827 + t599 * t665) * t675 - (t595 * t831 + t599 * t664) * t671) * t867 + (-(-t594 * t828 + t598 * t663) * t674 - (t594 * t832 + t598 * t662) * t670) * t874 + (-0.4e1 * t589 * t810 + t577) * t649 + (-0.4e1 * t588 * t812 + t576) * t648 + (-0.4e1 * t587 * t814 + t575) * t647 + (-0.4e1 * t586 * t822 + t574) * t646 + t892 * t637 + t893 * t636 + t894 * t635 + t895 * t634;];
taucX  = t1;
