% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V2G1A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:35:32
% EndTime: 2020-08-06 17:35:36
% DurationCPUTime: 4.13s
% Computational Cost: add. (22470->289), mult. (43956->530), div. (2541->10), fcn. (46854->22), ass. (0->228)
t769 = cos(qJ(2,3));
t778 = pkin(7) + pkin(6);
t729 = t769 * t778;
t763 = sin(qJ(2,3));
t708 = pkin(2) * t763 - t729;
t755 = sin(pkin(4));
t768 = cos(qJ(3,3));
t757 = cos(pkin(4));
t762 = sin(qJ(3,3));
t835 = t757 * t762;
t846 = t755 * t763;
t751 = t768 ^ 2;
t877 = pkin(3) * t751;
t663 = 0.1e1 / ((pkin(3) * t835 + t708 * t755) * t768 + pkin(2) * t835 + t846 * t877);
t771 = cos(qJ(2,2));
t730 = t771 * t778;
t765 = sin(qJ(2,2));
t709 = pkin(2) * t765 - t730;
t770 = cos(qJ(3,2));
t764 = sin(qJ(3,2));
t833 = t757 * t764;
t844 = t755 * t765;
t752 = t770 ^ 2;
t876 = pkin(3) * t752;
t664 = 0.1e1 / ((pkin(3) * t833 + t709 * t755) * t770 + pkin(2) * t833 + t844 * t876);
t773 = cos(qJ(2,1));
t731 = t773 * t778;
t767 = sin(qJ(2,1));
t710 = pkin(2) * t767 - t731;
t772 = cos(qJ(3,1));
t766 = sin(qJ(3,1));
t831 = t757 * t766;
t842 = t755 * t767;
t753 = t772 ^ 2;
t875 = pkin(3) * t753;
t665 = 0.1e1 / ((pkin(3) * t831 + t710 * t755) * t772 + pkin(2) * t831 + t842 * t875);
t758 = legFrame(3,3);
t734 = sin(t758);
t737 = cos(t758);
t754 = sin(pkin(8));
t756 = cos(pkin(8));
t684 = -t754 * t734 + t737 * t756;
t687 = t756 * t734 + t737 * t754;
t726 = t768 * pkin(3) + pkin(2);
t702 = t763 * t726 - t729;
t824 = t763 * t778;
t853 = (t726 * t769 + t824) * t757;
t657 = -t702 * t684 - t687 * t853;
t660 = -t684 * t853 + t702 * t687;
t705 = t726 * t835;
t841 = t755 * t768;
t669 = 0.1e1 / (t702 * t841 + t705);
t776 = xDP(2);
t777 = xDP(1);
t782 = 0.1e1 / pkin(3);
t642 = (t657 * t776 + t660 * t777) * t782 * t669;
t898 = -0.2e1 * t642;
t759 = legFrame(2,3);
t735 = sin(t759);
t738 = cos(t759);
t685 = -t754 * t735 + t738 * t756;
t688 = t756 * t735 + t738 * t754;
t727 = t770 * pkin(3) + pkin(2);
t703 = t765 * t727 - t730;
t822 = t765 * t778;
t852 = (t727 * t771 + t822) * t757;
t658 = -t703 * t685 - t688 * t852;
t661 = -t685 * t852 + t703 * t688;
t706 = t727 * t833;
t839 = t755 * t770;
t670 = 0.1e1 / (t703 * t839 + t706);
t643 = (t658 * t776 + t661 * t777) * t782 * t670;
t897 = -0.2e1 * t643;
t760 = legFrame(1,3);
t736 = sin(t760);
t739 = cos(t760);
t686 = -t754 * t736 + t739 * t756;
t689 = t756 * t736 + t739 * t754;
t728 = t772 * pkin(3) + pkin(2);
t704 = t767 * t728 - t731;
t820 = t767 * t778;
t851 = (t728 * t773 + t820) * t757;
t659 = -t704 * t686 - t689 * t851;
t662 = -t686 * t851 + t704 * t689;
t707 = t728 * t831;
t837 = t755 * t772;
t671 = 0.1e1 / (t704 * t837 + t707);
t644 = (t659 * t776 + t662 * t777) * t782 * t671;
t896 = -0.2e1 * t644;
t825 = t763 * t768;
t681 = pkin(3) * t825 + t708;
t834 = t757 * t763;
t651 = -t684 * t841 - (t684 * t834 + t769 * t687) * t762;
t654 = -t687 * t841 - (-t769 * t684 + t687 * t834) * t762;
t636 = (t651 * t777 + t654 * t776) * t663;
t864 = t636 * t778;
t799 = t762 * t864;
t829 = t757 * t782;
t847 = t755 * t762;
t783 = pkin(2) ^ 2;
t732 = t778 ^ 2 + t783;
t781 = pkin(3) ^ 2;
t880 = pkin(3) * t642;
t802 = t762 * t880;
t819 = 0.2e1 * pkin(2) * pkin(3);
t868 = (-t778 * t802 + (t751 * t781 + t768 * t819 + t732) * t636) * t636;
t606 = t663 * t829 * t868 + (-t757 * t799 + (-t681 * t847 + t757 * (pkin(2) * t768 + t877)) * t642) / (t681 * t841 + t705) * t642;
t627 = t799 - t880;
t615 = (-t768 * t868 - (pkin(2) * t642 - t627 * t768) * t880) * t663;
t633 = t636 ^ 2;
t639 = t642 ^ 2;
t871 = t762 * mrSges(3,1);
t717 = t768 * mrSges(3,2) + t871;
t796 = t768 * mrSges(3,1) - mrSges(3,2) * t762;
t672 = t717 * t846 - t757 * t796;
t733 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t745 = m(3) * pkin(2) + mrSges(2,1);
t749 = -m(1) - m(2) - m(3);
t865 = t636 * t769;
t892 = ((-t633 * t745 - t796 * (t633 + t639)) * t763 + (t636 * t733 + t717 * t898) * t865) * t755 + t672 * t606 + t749 * t615;
t823 = t765 * t770;
t683 = pkin(3) * t823 + t709;
t832 = t757 * t765;
t652 = -t685 * t839 - (t685 * t832 + t771 * t688) * t764;
t655 = -t688 * t839 - (-t771 * t685 + t688 * t832) * t764;
t637 = (t652 * t777 + t655 * t776) * t664;
t862 = t637 * t778;
t798 = t764 * t862;
t845 = t755 * t764;
t879 = pkin(3) * t643;
t801 = t764 * t879;
t867 = (-t778 * t801 + (t752 * t781 + t770 * t819 + t732) * t637) * t637;
t607 = t664 * t829 * t867 + (-t757 * t798 + (-t683 * t845 + t757 * (pkin(2) * t770 + t876)) * t643) / (t683 * t839 + t706) * t643;
t628 = t798 - t879;
t616 = (-t770 * t867 - (pkin(2) * t643 - t628 * t770) * t879) * t664;
t634 = t637 ^ 2;
t640 = t643 ^ 2;
t870 = t764 * mrSges(3,1);
t718 = t770 * mrSges(3,2) + t870;
t795 = t770 * mrSges(3,1) - mrSges(3,2) * t764;
t673 = t718 * t844 - t757 * t795;
t863 = t637 * t771;
t891 = ((-t634 * t745 - t795 * (t634 + t640)) * t765 + (t637 * t733 + t718 * t897) * t863) * t755 + t673 * t607 + t749 * t616;
t821 = t767 * t772;
t682 = pkin(3) * t821 + t710;
t830 = t757 * t767;
t653 = -t686 * t837 - (t686 * t830 + t773 * t689) * t766;
t656 = -t689 * t837 - (-t773 * t686 + t689 * t830) * t766;
t638 = (t653 * t777 + t656 * t776) * t665;
t860 = t638 * t778;
t797 = t766 * t860;
t843 = t755 * t766;
t878 = pkin(3) * t644;
t800 = t766 * t878;
t866 = (-t778 * t800 + (t753 * t781 + t772 * t819 + t732) * t638) * t638;
t608 = t665 * t829 * t866 + (-t757 * t797 + (-t682 * t843 + t757 * (pkin(2) * t772 + t875)) * t644) / (t682 * t837 + t707) * t644;
t629 = t797 - t878;
t617 = (-t772 * t866 - (pkin(2) * t644 - t629 * t772) * t878) * t665;
t635 = t638 ^ 2;
t641 = t644 ^ 2;
t869 = t766 * mrSges(3,1);
t719 = t772 * mrSges(3,2) + t869;
t794 = t772 * mrSges(3,1) - mrSges(3,2) * t766;
t674 = t719 * t842 - t757 * t794;
t861 = t638 * t773;
t890 = ((-t635 * t745 - t794 * (t635 + t641)) * t767 + (t638 * t733 + t719 * t896) * t861) * t755 + t674 * t608 + t749 * t617;
t713 = pkin(2) * t773 + t820;
t785 = pkin(3) * t843 - t710 * t757;
t889 = t713 * t756 + t754 * t785;
t712 = pkin(2) * t771 + t822;
t786 = pkin(3) * t845 - t709 * t757;
t888 = t712 * t756 + t754 * t786;
t711 = pkin(2) * t769 + t824;
t787 = pkin(3) * t847 - t708 * t757;
t887 = t711 * t756 + t754 * t787;
t761 = Ifges(3,1) - Ifges(3,2);
t774 = mrSges(3,2) * pkin(2);
t883 = -2 * Ifges(3,4);
t886 = Ifges(3,4) + t753 * t883 + (-t761 * t766 + t774) * t772;
t885 = Ifges(3,4) + t752 * t883 + (-t761 * t764 + t774) * t770;
t884 = Ifges(3,4) + t751 * t883 + (-t761 * t762 + t774) * t768;
t775 = mrSges(3,1) * pkin(2);
t740 = pkin(6) * mrSges(3,2) - Ifges(3,6);
t882 = -t740 / 0.2e1;
t741 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t881 = t741 / 0.2e1;
t859 = t639 * t717;
t858 = t640 * t718;
t857 = t641 * t719;
t675 = (t796 + t745) * t769 + t763 * t733;
t856 = t675 * t755;
t676 = (t795 + t745) * t771 + t765 * t733;
t855 = t676 * t755;
t677 = (t794 + t745) * t773 + t767 * t733;
t854 = t677 * t755;
t840 = t755 * t769;
t838 = t755 * t771;
t836 = t755 * t773;
t600 = (((t636 * t840 + t757 * t642) * t877 + ((-t802 + t864) * t763 + pkin(2) * t865) * t841 + t627 * t757) * t636 + (t642 * t840 + (t751 * t757 - t825 * t847 - t757) * t636) * t880) * t663;
t818 = -t600 * t856 - t757 * t859 + t892;
t601 = (((t637 * t838 + t757 * t643) * t876 + ((-t801 + t862) * t765 + pkin(2) * t863) * t839 + t628 * t757) * t637 + (t643 * t838 + (t752 * t757 - t823 * t845 - t757) * t637) * t879) * t664;
t817 = -t601 * t855 - t757 * t858 + t891;
t602 = (((t638 * t836 + t757 * t644) * t875 + ((-t800 + t860) * t767 + pkin(2) * t861) * t837 + t629 * t757) * t638 + (t644 * t836 + (t753 * t757 - t821 * t843 - t757) * t638) * t878) * t665;
t816 = -t602 * t854 - t757 * t857 + t890;
t809 = 0.2e1 * t774;
t805 = pkin(2) * t847;
t804 = pkin(2) * t845;
t803 = pkin(2) * t843;
t696 = t740 * t768 + t762 * t741;
t784 = -0.2e1 * pkin(6) * mrSges(3,3) + (-pkin(6) ^ 2 - t783) * m(3) - Ifges(3,1) - Ifges(2,3);
t793 = ((t762 * t882 + t768 * t881) * t642 + (t775 * t762 + t884) * t636) * t898 - t615 * t856 + (t761 * t751 - 0.2e1 * (Ifges(3,4) * t762 + t775) * t768 + t762 * t809 + t784) * t600 + t696 * t606;
t697 = t740 * t770 + t764 * t741;
t792 = ((t764 * t882 + t770 * t881) * t643 + (t775 * t764 + t885) * t637) * t897 - t616 * t855 + (t761 * t752 - 0.2e1 * (Ifges(3,4) * t764 + t775) * t770 + t764 * t809 + t784) * t601 + t697 * t607;
t698 = t740 * t772 + t766 * t741;
t791 = ((t766 * t882 + t772 * t881) * t644 + (t775 * t766 + t886) * t638) * t896 - t617 * t854 + (t761 * t753 - 0.2e1 * (Ifges(3,4) * t766 + t775) * t772 + t766 * t809 + t784) * t602 + t698 * t608;
t790 = (-Ifges(3,3) * t606 + t696 * t600 + t672 * t615 + t633 * (pkin(2) * t871 + t884)) * t669;
t789 = (-Ifges(3,3) * t607 + t697 * t601 + t673 * t616 + t634 * (pkin(2) * t870 + t885)) * t670;
t788 = (-Ifges(3,3) * t608 + t698 * t602 + t674 * t617 + t635 * (pkin(2) * t869 + t886)) * t671;
t695 = t754 * t773 + t756 * t830;
t694 = t754 * t771 + t756 * t832;
t693 = t754 * t769 + t756 * t834;
t692 = t754 * t830 - t756 * t773;
t691 = t754 * t832 - t756 * t771;
t690 = t754 * t834 - t756 * t769;
t668 = t754 * t713 - t756 * t785;
t667 = t754 * t712 - t756 * t786;
t666 = t754 * t711 - t756 * t787;
t1 = [(t791 * t653 + t816 * (-(t692 * t739 + t736 * t695) * t875 + (-t668 * t736 + t889 * t739) * t772 + t689 * t803)) * t665 + (t792 * t652 + t817 * (-(t691 * t738 + t735 * t694) * t876 + (-t667 * t735 + t888 * t738) * t770 + t688 * t804)) * t664 + (t793 * t651 + t818 * (-(t690 * t737 + t734 * t693) * t877 + (-t666 * t734 + t887 * t737) * t768 + t687 * t805)) * t663 + (t660 * t790 + t661 * t789 + t662 * t788) * t782; (t791 * t656 + t816 * ((-t736 * t692 + t695 * t739) * t875 + (t668 * t739 + t889 * t736) * t772 - t686 * t803)) * t665 + (t792 * t655 + t817 * ((-t735 * t691 + t694 * t738) * t876 + (t667 * t738 + t888 * t735) * t770 - t685 * t804)) * t664 + (t793 * t654 + t818 * ((-t734 * t690 + t693 * t737) * t877 + (t666 * t737 + t887 * t734) * t768 - t684 * t805)) * t663 + (t657 * t790 + t658 * t789 + t659 * t788) * t782; (-t857 - t858 - t859) * t757 + (-t600 * t675 - t601 * t676 - t602 * t677) * t755 + t890 + t891 + t892;];
taucX  = t1;
