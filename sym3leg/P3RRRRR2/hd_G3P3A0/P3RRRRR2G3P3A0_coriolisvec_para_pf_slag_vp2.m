% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRRRR2G3P3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 21:13
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRRRR2G3P3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G3P3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G3P3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G3P3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G3P3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G3P3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR2G3P3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR2G3P3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G3P3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G3P3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:11:58
% EndTime: 2020-03-09 21:12:01
% DurationCPUTime: 3.11s
% Computational Cost: add. (9273->249), mult. (23730->506), div. (9792->11), fcn. (27528->27), ass. (0->239)
t756 = legFrame(1,2);
t736 = sin(t756);
t739 = cos(t756);
t772 = cos(qJ(3,1));
t763 = sin(qJ(3,1));
t773 = cos(qJ(2,1));
t892 = t773 * pkin(1);
t829 = t763 * t892;
t764 = sin(qJ(2,1));
t765 = sin(qJ(1,1));
t774 = cos(qJ(1,1));
t718 = t765 * t764 - t774 * t773;
t749 = t772 ^ 2;
t832 = pkin(2) * t718 * t749;
t864 = t739 * t763;
t895 = pkin(1) * t774;
t679 = -t736 * t832 + (-pkin(2) * t864 + t736 * t895) * t772 - t739 * t829;
t776 = xDP(2);
t751 = 0.1e1 / t772 ^ 2;
t742 = 0.1e1 / t764;
t780 = 0.1e1 / pkin(2);
t782 = 0.1e1 / pkin(1);
t841 = t780 * t782;
t811 = t742 * t841;
t799 = t751 * t811;
t673 = t679 * t776 * t799;
t867 = t736 * t763;
t682 = t739 * t832 + (-pkin(2) * t867 - t739 * t895) * t772 - t736 * t829;
t777 = xDP(1);
t676 = t682 * t777 * t799;
t712 = pkin(2) * (t774 * t764 + t765 * t773) * t772 + t765 * pkin(1);
t775 = xDP(3);
t750 = 0.1e1 / t772;
t800 = t750 * t811;
t694 = t712 * t775 * t800;
t664 = t676 + t673 + t694;
t873 = t718 * t772;
t699 = t736 * t873 + t864;
t700 = -t739 * t873 + t867;
t861 = t742 * t782;
t812 = t750 * t861;
t842 = t775 * t782;
t870 = sin(qJ(1,1) + qJ(2,1)) * t742;
t670 = -t842 * t870 + (t699 * t776 + t700 * t777) * t812;
t835 = -t664 - t670;
t658 = t835 ^ 2;
t755 = legFrame(2,2);
t735 = sin(t755);
t738 = cos(t755);
t769 = cos(qJ(3,2));
t760 = sin(qJ(3,2));
t770 = cos(qJ(2,2));
t893 = t770 * pkin(1);
t830 = t760 * t893;
t761 = sin(qJ(2,2));
t762 = sin(qJ(1,2));
t771 = cos(qJ(1,2));
t717 = t762 * t761 - t771 * t770;
t746 = t769 ^ 2;
t833 = pkin(2) * t717 * t746;
t865 = t738 * t760;
t896 = pkin(1) * t771;
t678 = -t735 * t833 + (-pkin(2) * t865 + t735 * t896) * t769 - t738 * t830;
t748 = 0.1e1 / t769 ^ 2;
t741 = 0.1e1 / t761;
t813 = t741 * t841;
t801 = t748 * t813;
t672 = t678 * t776 * t801;
t868 = t735 * t760;
t681 = t738 * t833 + (-pkin(2) * t868 - t738 * t896) * t769 - t735 * t830;
t675 = t681 * t777 * t801;
t711 = pkin(2) * (t771 * t761 + t762 * t770) * t769 + t762 * pkin(1);
t747 = 0.1e1 / t769;
t802 = t747 * t813;
t693 = t711 * t775 * t802;
t663 = t675 + t672 + t693;
t874 = t717 * t769;
t697 = t735 * t874 + t865;
t698 = -t738 * t874 + t868;
t862 = t741 * t782;
t814 = t747 * t862;
t871 = sin(qJ(1,2) + qJ(2,2)) * t741;
t669 = -t842 * t871 + (t697 * t776 + t698 * t777) * t814;
t836 = -t663 - t669;
t657 = t836 ^ 2;
t754 = legFrame(3,2);
t734 = sin(t754);
t737 = cos(t754);
t766 = cos(qJ(3,3));
t757 = sin(qJ(3,3));
t767 = cos(qJ(2,3));
t894 = t767 * pkin(1);
t831 = t757 * t894;
t758 = sin(qJ(2,3));
t759 = sin(qJ(1,3));
t768 = cos(qJ(1,3));
t716 = t759 * t758 - t768 * t767;
t743 = t766 ^ 2;
t834 = pkin(2) * t716 * t743;
t866 = t737 * t757;
t897 = pkin(1) * t768;
t677 = -t734 * t834 + (-pkin(2) * t866 + t734 * t897) * t766 - t737 * t831;
t745 = 0.1e1 / t766 ^ 2;
t740 = 0.1e1 / t758;
t815 = t740 * t841;
t803 = t745 * t815;
t671 = t677 * t776 * t803;
t869 = t734 * t757;
t680 = t737 * t834 + (-pkin(2) * t869 - t737 * t897) * t766 - t734 * t831;
t674 = t680 * t777 * t803;
t710 = pkin(2) * (t768 * t758 + t759 * t767) * t766 + t759 * pkin(1);
t744 = 0.1e1 / t766;
t804 = t744 * t815;
t692 = t710 * t775 * t804;
t662 = t674 + t671 + t692;
t875 = t716 * t766;
t695 = t734 * t875 + t866;
t696 = -t737 * t875 + t869;
t863 = t740 * t782;
t816 = t744 * t863;
t872 = sin(qJ(1,3) + qJ(2,3)) * t740;
t668 = -t842 * t872 + (t695 * t776 + t696 * t777) * t816;
t837 = -t662 - t668;
t656 = t837 ^ 2;
t907 = 0.2e1 * pkin(1);
t665 = t668 ^ 2;
t906 = t665 * t894;
t666 = t669 ^ 2;
t905 = t666 * t893;
t667 = t670 ^ 2;
t904 = t667 * t892;
t903 = -0.2e1 * t766;
t902 = -0.2e1 * t769;
t901 = -0.2e1 * t772;
t900 = pkin(1) * t758;
t899 = pkin(1) * t761;
t898 = pkin(1) * t764;
t891 = -Ifges(3,1) - Ifges(2,3);
t890 = Ifges(3,4) * t757;
t889 = Ifges(3,4) * t760;
t888 = Ifges(3,4) * t763;
t707 = (t734 * t777 + t737 * t776) * t780 * t744;
t704 = t707 ^ 2;
t887 = (t704 / 0.2e1 + (t668 + t662 / 0.2e1) * t662) * t758;
t708 = (t735 * t777 + t738 * t776) * t780 * t747;
t705 = t708 ^ 2;
t886 = (t705 / 0.2e1 + (t669 + t663 / 0.2e1) * t663) * t761;
t709 = (t736 * t777 + t739 * t776) * t780 * t750;
t706 = t709 ^ 2;
t885 = (t706 / 0.2e1 + (t670 + t664 / 0.2e1) * t664) * t764;
t884 = t837 * t707;
t883 = t837 * t743;
t882 = t836 * t708;
t881 = t836 * t746;
t880 = t835 * t709;
t879 = t835 * t749;
t878 = t704 * t757;
t877 = t705 * t760;
t876 = t706 * t763;
t752 = Ifges(3,2) - Ifges(3,1);
t860 = t752 * t743;
t859 = t752 * t746;
t858 = t752 * t749;
t857 = t752 * t757;
t856 = t752 * t760;
t855 = t752 * t763;
t753 = mrSges(2,2) - mrSges(3,3);
t854 = t753 * t767;
t853 = t753 * t770;
t852 = t753 * t773;
t851 = t757 * t758;
t850 = t760 * t761;
t849 = t763 * t764;
t848 = t766 * t767;
t847 = t767 * t707;
t846 = t769 * t770;
t845 = t770 * t708;
t844 = t772 * t773;
t843 = t773 * t709;
t653 = t674 / 0.2e1 + t671 / 0.2e1 + t692 / 0.2e1 + t668;
t779 = pkin(2) ^ 2;
t819 = t707 * t851;
t629 = (-t779 * t883 + (pkin(1) * t668 + (0.2e1 * t653 * t848 - t819) * pkin(2)) * pkin(1)) * t668 * t804 + (-pkin(2) * t883 + (-t837 * t848 - t819) * pkin(1)) * t662 * t816 + ((pkin(1) * t837 * t851 + pkin(2) * t707) * t766 + pkin(1) * t847) * t745 * t707 * t863;
t638 = (-t906 + (-t766 * t656 - t704 * t744) * pkin(2)) * t863;
t728 = mrSges(3,1) * t894;
t722 = t757 * mrSges(3,2) - mrSges(2,1);
t788 = (t722 * t767 + t753 * t758) * pkin(1);
t798 = -t860 + t891;
t683 = -(t728 + 0.2e1 * t890) * t766 + t788 + t798;
t713 = -t766 * (-mrSges(3,2) * t900 + Ifges(3,6)) + t757 * (mrSges(3,1) * t900 - Ifges(3,5));
t785 = (t704 * Ifges(3,5) + 0.2e1 * t857 * t884) * t766 - Ifges(3,6) * t878 - (0.4e1 * t743 - 0.2e1) * Ifges(3,4) * t884;
t789 = -(m(2) + m(3)) * pkin(1) ^ 2 - Ifges(1,3) + t891;
t822 = t744 * t878;
t825 = t837 * t847;
t840 = t785 + ((-mrSges(3,1) * t887 + mrSges(3,2) * t825) * t766 + (mrSges(3,1) * t825 + mrSges(3,2) * t887) * t757 + (-mrSges(2,1) * t758 - t854) * t662 * t653) * t907 + (-t860 + (t728 + t890) * t903 + 0.2e1 * t788 + t789) * t638 + t683 * t629 - t713 * t822;
t654 = t675 / 0.2e1 + t672 / 0.2e1 + t693 / 0.2e1 + t669;
t818 = t708 * t850;
t630 = (-t779 * t881 + (pkin(1) * t669 + (0.2e1 * t654 * t846 - t818) * pkin(2)) * pkin(1)) * t669 * t802 + (-pkin(2) * t881 + (-t836 * t846 - t818) * pkin(1)) * t663 * t814 + ((pkin(1) * t836 * t850 + pkin(2) * t708) * t769 + pkin(1) * t845) * t748 * t708 * t862;
t639 = (-t905 + (-t769 * t657 - t705 * t747) * pkin(2)) * t862;
t729 = mrSges(3,1) * t893;
t723 = t760 * mrSges(3,2) - mrSges(2,1);
t787 = (t723 * t770 + t753 * t761) * pkin(1);
t797 = -t859 + t891;
t684 = -(t729 + 0.2e1 * t889) * t769 + t787 + t797;
t714 = -t769 * (-mrSges(3,2) * t899 + Ifges(3,6)) + t760 * (mrSges(3,1) * t899 - Ifges(3,5));
t784 = (t705 * Ifges(3,5) + 0.2e1 * t856 * t882) * t769 - Ifges(3,6) * t877 - (0.4e1 * t746 - 0.2e1) * Ifges(3,4) * t882;
t821 = t747 * t877;
t824 = t836 * t845;
t839 = t784 + ((-mrSges(3,1) * t886 + mrSges(3,2) * t824) * t769 + (mrSges(3,1) * t824 + mrSges(3,2) * t886) * t760 + (-mrSges(2,1) * t761 - t853) * t663 * t654) * t907 + (-t859 + (t729 + t889) * t902 + 0.2e1 * t787 + t789) * t639 + t684 * t630 - t714 * t821;
t655 = t676 / 0.2e1 + t673 / 0.2e1 + t694 / 0.2e1 + t670;
t817 = t709 * t849;
t631 = (-t779 * t879 + (pkin(1) * t670 + (0.2e1 * t655 * t844 - t817) * pkin(2)) * pkin(1)) * t670 * t800 + (-pkin(2) * t879 + (-t835 * t844 - t817) * pkin(1)) * t664 * t812 + ((pkin(1) * t835 * t849 + pkin(2) * t709) * t772 + pkin(1) * t843) * t751 * t709 * t861;
t640 = (-t904 + (-t772 * t658 - t706 * t750) * pkin(2)) * t861;
t730 = mrSges(3,1) * t892;
t724 = t763 * mrSges(3,2) - mrSges(2,1);
t786 = (t724 * t773 + t753 * t764) * pkin(1);
t796 = -t858 + t891;
t685 = -(t730 + 0.2e1 * t888) * t772 + t786 + t796;
t715 = -t772 * (-mrSges(3,2) * t898 + Ifges(3,6)) + t763 * (mrSges(3,1) * t898 - Ifges(3,5));
t783 = (t706 * Ifges(3,5) + 0.2e1 * t855 * t880) * t772 - Ifges(3,6) * t876 - (0.4e1 * t749 - 0.2e1) * Ifges(3,4) * t880;
t820 = t750 * t876;
t823 = t835 * t843;
t838 = t783 + ((-mrSges(3,1) * t885 + mrSges(3,2) * t823) * t772 + (mrSges(3,1) * t823 + mrSges(3,2) * t885) * t763 + (-mrSges(2,1) * t764 - t852) * t664 * t655) * t907 + (-t858 + (t730 + t888) * t901 + 0.2e1 * t786 + t789) * t640 + t685 * t631 - t715 * t820;
t719 = -Ifges(3,5) * t757 - Ifges(3,6) * t766;
t810 = t740 * (t683 * t638 + (t890 * t903 + t798) * t629 - t719 * t822 + (t854 + (mrSges(3,1) * t766 - t722) * t758) * t665 * pkin(1) + t785);
t720 = -Ifges(3,5) * t760 - Ifges(3,6) * t769;
t809 = t741 * (t684 * t639 + (t889 * t902 + t797) * t630 - t720 * t821 + (t853 + (mrSges(3,1) * t769 - t723) * t761) * t666 * pkin(1) + t784);
t721 = -Ifges(3,5) * t763 - Ifges(3,6) * t772;
t808 = t742 * (t685 * t640 + (t888 * t901 + t796) * t631 - t721 * t820 + (t852 + (mrSges(3,1) * t772 - t724) * t764) * t667 * pkin(1) + t783);
t807 = (Ifges(3,3) * t822 + t719 * t629 + t713 * t638 + (mrSges(3,2) * t906 + t656 * t857) * t766 + mrSges(3,1) * t665 * t831 + (-0.2e1 * t743 + 0.1e1) * t656 * Ifges(3,4)) * t744;
t806 = (Ifges(3,3) * t821 + t720 * t630 + t714 * t639 + (mrSges(3,2) * t905 + t657 * t856) * t769 + mrSges(3,1) * t666 * t830 + (-0.2e1 * t746 + 0.1e1) * t657 * Ifges(3,4)) * t747;
t805 = (Ifges(3,3) * t820 + t721 * t631 + t715 * t640 + (mrSges(3,2) * t904 + t658 * t855) * t772 + mrSges(3,1) * t667 * t829 + (-0.2e1 * t749 + 0.1e1) * t658 * Ifges(3,4)) * t750;
t795 = t840 * t744 * t740;
t794 = t839 * t747 * t741;
t793 = t838 * t750 * t742;
t792 = t745 * t810;
t791 = t748 * t809;
t790 = t751 * t808;
t1 = [(t696 * t795 + t698 * t794 + t700 * t793) * t782 + (t736 * t805 + t735 * t806 + t734 * t807 + (t680 * t792 + t681 * t791 + t682 * t790) * t782) * t780; (t695 * t795 + t697 * t794 + t699 * t793) * t782 + (t739 * t805 + t738 * t806 + t737 * t807 + (t677 * t792 + t678 * t791 + t679 * t790) * t782) * t780; (-t838 * t870 - t839 * t871 - t840 * t872 + (t744 * t710 * t810 + t747 * t711 * t809 + t750 * t712 * t808) * t780) * t782;];
taucX  = t1;
