% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR12V2G1A0
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2020-08-06 19:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:15:05
% EndTime: 2020-08-06 19:15:13
% DurationCPUTime: 8.07s
% Computational Cost: add. (60408->446), mult. (73062->691), div. (8943->6), fcn. (64485->18), ass. (0->294)
t821 = (pkin(2) + pkin(3));
t986 = -2 * t821;
t820 = (pkin(5) - pkin(6));
t990 = t820 * t821;
t878 = pkin(2) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t981 = m(3) * pkin(2) + mrSges(3,1);
t989 = qJ(3,1) * t981 + t878;
t988 = qJ(3,2) * t981 + t878;
t987 = qJ(3,3) * t981 + t878;
t985 = (-pkin(5) / 0.2e1 + pkin(6) / 0.2e1) * t986;
t827 = qJ(3,1) ^ 2;
t830 = pkin(2) ^ 2;
t984 = (t827 - t830) * m(3);
t825 = qJ(3,2) ^ 2;
t983 = (t825 - t830) * m(3);
t823 = qJ(3,3) ^ 2;
t982 = (t823 - t830) * m(3);
t762 = qJ(3,1) * m(3) + mrSges(3,3);
t760 = m(3) * qJ(3,3) + mrSges(3,3);
t761 = m(3) * qJ(3,2) + mrSges(3,3);
t796 = legFrame(3,3);
t763 = sin(t796);
t766 = cos(t796);
t804 = sin(qJ(1,3));
t810 = cos(qJ(1,3));
t708 = t763 * t804 - t766 * t810;
t709 = t763 * t810 + t766 * t804;
t803 = sin(qJ(2,3));
t913 = t803 * qJ(3,3);
t742 = pkin(1) + t913;
t809 = cos(qJ(2,3));
t903 = t821 * t809;
t723 = 0.1e1 / (t742 + t903);
t818 = xDP(2);
t819 = xDP(1);
t698 = (-t708 * t818 - t709 * t819) * t723;
t789 = t809 ^ 2;
t980 = t698 * t789;
t797 = legFrame(2,3);
t764 = sin(t797);
t767 = cos(t797);
t806 = sin(qJ(1,2));
t812 = cos(qJ(1,2));
t710 = t764 * t806 - t767 * t812;
t711 = t764 * t812 + t767 * t806;
t805 = sin(qJ(2,2));
t911 = t805 * qJ(3,2);
t744 = pkin(1) + t911;
t811 = cos(qJ(2,2));
t902 = t821 * t811;
t724 = 0.1e1 / (t744 + t902);
t699 = (-t710 * t818 - t711 * t819) * t724;
t791 = t811 ^ 2;
t979 = t699 * t791;
t798 = legFrame(1,3);
t765 = sin(t798);
t768 = cos(t798);
t808 = sin(qJ(1,1));
t814 = cos(qJ(1,1));
t712 = t765 * t808 - t768 * t814;
t713 = t765 * t814 + t768 * t808;
t807 = sin(qJ(2,1));
t909 = t807 * qJ(3,1);
t746 = pkin(1) + t909;
t813 = cos(qJ(2,1));
t901 = t821 * t813;
t725 = 0.1e1 / (t746 + t901);
t700 = (-t712 * t818 - t713 * t819) * t725;
t793 = t813 ^ 2;
t978 = t700 * t793;
t950 = pkin(1) * t803;
t747 = qJ(3,3) + t950;
t753 = t820 * t810;
t704 = t804 * t747 - t803 * t753;
t741 = pkin(1) + 0.2e1 * t913;
t714 = t804 * t741 - t753;
t750 = t804 * t820;
t841 = t747 * t810 + t803 * t750;
t855 = t741 * t810 + t750;
t916 = (qJ(3,3) + t821) * (-qJ(3,3) + t821);
t859 = t789 * t916;
t683 = -t708 * t859 - (t763 * t714 - t855 * t766) * t903 - (t763 * t704 - t841 * t766) * qJ(3,3);
t684 = t709 * t859 + (t714 * t766 + t855 * t763) * t903 + (t704 * t766 + t841 * t763) * qJ(3,3);
t906 = t821 * t803;
t733 = -t809 * qJ(3,3) + t906;
t817 = xDP(3);
t824 = 0.1e1 / qJ(3,3);
t923 = t723 * t824;
t662 = t733 * t824 * t817 + (t683 * t819 + t684 * t818) * t923;
t715 = t804 * t742 - t753;
t854 = t742 * t810 + t750;
t689 = t708 * t903 + t715 * t763 - t854 * t766;
t690 = t709 * t903 + t715 * t766 + t763 * t854;
t674 = (t803 * t817 + (-t689 * t819 + t690 * t818) * t809 * t723) * t824;
t974 = t821 * t674;
t648 = t662 - t974;
t641 = t821 * t648 - t823 * t674;
t656 = t674 * t820;
t659 = t662 * t820;
t831 = pkin(1) ^ 2;
t887 = -pkin(5) ^ 2 - t831;
t967 = -0.2e1 * pkin(5);
t736 = (t967 + pkin(6)) * pkin(6) - t887;
t788 = t809 * t789;
t931 = t698 * t820;
t864 = t803 * t931;
t865 = t698 * t923;
t868 = t674 * t923;
t871 = t662 * t923;
t879 = pkin(1) * qJ(3,3) * t698;
t886 = 0.2e1 * pkin(1);
t920 = t736 * t821;
t930 = t698 * t821;
t932 = t698 * t803;
t938 = (t864 - t974) * t809;
t956 = t789 - 0.1e1;
t962 = -0.3e1 * t823;
t965 = -0.2e1 * t789;
t966 = t821 ^ 2;
t971 = t823 - t966;
t629 = (-(t966 + t962) * t788 * t930 + (-0.3e1 * (-t823 / 0.3e1 + t966) * t913 + t971 * t886) * t980 + (((-t659 - 0.4e1 * t879) * t821 - t971 * t656) * t803 + t930 * t962 - t698 * t920) * t809) * t865 + ((t641 * t821 + t864 * t916) * t809 + t641 * pkin(1)) * t868 + (pkin(1) * t674 - t938) * t821 * t871 + (((t656 * t986 + t659) * t789 + (-t736 - t823) * t932 - 0.2e1 * t879 - t659 + t674 * t985) * t865 + (t641 * t803 + (t965 + 0.1e1) * t698 * t990) * t868 + (t674 * t906 + t956 * t931) * t871) * qJ(3,3);
t872 = qJ(3,3) * t789 * t820;
t882 = 2 * t821;
t885 = -0.2e1 * t950;
t891 = -(t820 ^ 2) - t831;
t912 = t803 * t809;
t632 = (-t674 * t872 + (t674 * t990 - t659) * t912 + (-t788 * t916 + (qJ(3,3) * t885 - t823 + t891) * t809 - t742 * t789 * t882) * t698) * t865 + (-t698 * t872 + (t648 + t864) * t903 + t648 * t742) * t868 + (t674 * t742 - t938) * t871;
t638 = (0.2e1 * t662 * t803 - 0.2e1 * t733 * t674 + t931) * t698 * t723;
t695 = t698 ^ 2;
t759 = mrSges(2,1) + t981;
t726 = pkin(2) * mrSges(3,2) + t759 * pkin(5) - Ifges(3,4) - Ifges(2,5);
t756 = -mrSges(2,2) + t760;
t945 = Ifges(3,6) - Ifges(2,6);
t847 = -qJ(3,3) * mrSges(3,2) + t945;
t701 = -t809 * (t756 * pkin(5) - t847) + t726 * t803;
t942 = mrSges(3,3) * qJ(3,3);
t780 = 0.2e1 * t942;
t781 = -0.2e1 * t942;
t957 = mrSges(3,1) * pkin(2);
t787 = -0.2e1 * t957;
t946 = Ifges(2,1) + Ifges(3,1);
t856 = Ifges(2,2) + Ifges(3,3) - t946;
t846 = t787 - t856;
t877 = t787 - Ifges(3,2) - Ifges(2,3);
t941 = t662 * t760;
t947 = t759 * pkin(1);
t953 = pkin(1) * t756;
t897 = t701 * t638 + (-m(3) * (t823 + t830) + t781 + t877) * t632 + t981 * t629 + 0.2e1 * t674 * t941 + (t987 * t965 - ((t780 + t846 + t982) * t803 + t953) * t809 + t803 * t947 + t987) * t695;
t977 = t809 * t897;
t949 = pkin(1) * t805;
t748 = qJ(3,2) + t949;
t754 = t820 * t812;
t705 = t806 * t748 - t805 * t754;
t743 = pkin(1) + 0.2e1 * t911;
t716 = t806 * t743 - t754;
t751 = t806 * t820;
t840 = t748 * t812 + t805 * t751;
t853 = t743 * t812 + t751;
t915 = (qJ(3,2) + t821) * (-qJ(3,2) + t821);
t858 = t791 * t915;
t685 = -t710 * t858 - (t764 * t716 - t853 * t767) * t902 - (t764 * t705 - t840 * t767) * qJ(3,2);
t686 = t711 * t858 + (t716 * t767 + t853 * t764) * t902 + (t705 * t767 + t840 * t764) * qJ(3,2);
t905 = t821 * t805;
t734 = -t811 * qJ(3,2) + t905;
t826 = 0.1e1 / qJ(3,2);
t922 = t724 * t826;
t663 = t734 * t826 * t817 + (t685 * t819 + t686 * t818) * t922;
t717 = t806 * t744 - t754;
t852 = t744 * t812 + t751;
t691 = t710 * t902 + t717 * t764 - t852 * t767;
t692 = t711 * t902 + t717 * t767 + t764 * t852;
t675 = (t805 * t817 + (-t691 * t819 + t692 * t818) * t811 * t724) * t826;
t973 = t821 * t675;
t649 = t663 - t973;
t642 = t821 * t649 - t825 * t675;
t657 = t675 * t820;
t660 = t663 * t820;
t790 = t811 * t791;
t928 = t699 * t820;
t862 = t805 * t928;
t863 = t699 * t922;
t867 = t675 * t922;
t870 = t663 * t922;
t880 = pkin(1) * qJ(3,2) * t699;
t927 = t699 * t821;
t929 = t699 * t805;
t937 = (t862 - t973) * t811;
t955 = t791 - 0.1e1;
t961 = -0.3e1 * t825;
t964 = -0.2e1 * t791;
t970 = t825 - t966;
t630 = (-(t966 + t961) * t790 * t927 + (-0.3e1 * (-t825 / 0.3e1 + t966) * t911 + t970 * t886) * t979 + (((-t660 - 0.4e1 * t880) * t821 - t970 * t657) * t805 + t927 * t961 - t699 * t920) * t811) * t863 + ((t642 * t821 + t862 * t915) * t811 + t642 * pkin(1)) * t867 + (pkin(1) * t675 - t937) * t821 * t870 + (((t657 * t986 + t660) * t791 + (-t736 - t825) * t929 - 0.2e1 * t880 - t660 + t675 * t985) * t863 + (t642 * t805 + (t964 + 0.1e1) * t699 * t990) * t867 + (t675 * t905 + t955 * t928) * t870) * qJ(3,2);
t873 = qJ(3,2) * t791 * t820;
t884 = -0.2e1 * t949;
t910 = t805 * t811;
t633 = (-t675 * t873 + (t675 * t990 - t660) * t910 + (-t790 * t915 + (qJ(3,2) * t884 - t825 + t891) * t811 - t744 * t791 * t882) * t699) * t863 + (-t699 * t873 + (t649 + t862) * t902 + t649 * t744) * t867 + (t675 * t744 - t937) * t870;
t639 = (0.2e1 * t805 * t663 - 0.2e1 * t734 * t675 + t928) * t699 * t724;
t696 = t699 ^ 2;
t757 = -mrSges(2,2) + t761;
t848 = -qJ(3,2) * mrSges(3,2) + t945;
t702 = -t811 * (t757 * pkin(5) - t848) + t726 * t805;
t943 = mrSges(3,3) * qJ(3,2);
t782 = 0.2e1 * t943;
t783 = -0.2e1 * t943;
t940 = t663 * t761;
t952 = pkin(1) * t757;
t896 = t702 * t639 + (-m(3) * (t825 + t830) + t783 + t877) * t633 + t981 * t630 + 0.2e1 * t675 * t940 + (t988 * t964 - ((t782 + t846 + t983) * t805 + t952) * t811 + t805 * t947 + t988) * t696;
t976 = t811 * t896;
t948 = pkin(1) * t807;
t749 = qJ(3,1) + t948;
t755 = t820 * t814;
t706 = t808 * t749 - t807 * t755;
t745 = pkin(1) + 0.2e1 * t909;
t718 = t808 * t745 - t755;
t752 = t808 * t820;
t839 = t749 * t814 + t807 * t752;
t851 = t745 * t814 + t752;
t914 = (qJ(3,1) + t821) * (-qJ(3,1) + t821);
t857 = t793 * t914;
t687 = -t712 * t857 - (t765 * t718 - t851 * t768) * t901 - (t765 * t706 - t839 * t768) * qJ(3,1);
t688 = t713 * t857 + (t718 * t768 + t851 * t765) * t901 + (t706 * t768 + t839 * t765) * qJ(3,1);
t904 = t821 * t807;
t735 = -t813 * qJ(3,1) + t904;
t828 = 0.1e1 / qJ(3,1);
t921 = t725 * t828;
t664 = t735 * t828 * t817 + (t687 * t819 + t688 * t818) * t921;
t719 = t808 * t746 - t755;
t850 = t746 * t814 + t752;
t693 = t712 * t901 + t719 * t765 - t850 * t768;
t694 = t713 * t901 + t719 * t768 + t765 * t850;
t676 = (t807 * t817 + (-t693 * t819 + t694 * t818) * t813 * t725) * t828;
t972 = t821 * t676;
t647 = t664 - t972;
t643 = t821 * t647 - t827 * t676;
t658 = t676 * t820;
t661 = t664 * t820;
t792 = t813 * t793;
t925 = t700 * t820;
t860 = t807 * t925;
t861 = t700 * t921;
t866 = t676 * t921;
t869 = t664 * t921;
t881 = pkin(1) * qJ(3,1) * t700;
t924 = t700 * t821;
t926 = t700 * t807;
t936 = (t860 - t972) * t813;
t954 = t793 - 0.1e1;
t960 = -0.3e1 * t827;
t963 = -0.2e1 * t793;
t969 = t827 - t966;
t631 = (-(t966 + t960) * t792 * t924 + (-0.3e1 * (-t827 / 0.3e1 + t966) * t909 + t969 * t886) * t978 + (((-t661 - 0.4e1 * t881) * t821 - t969 * t658) * t807 + t924 * t960 - t700 * t920) * t813) * t861 + ((t643 * t821 + t860 * t914) * t813 + t643 * pkin(1)) * t866 + (pkin(1) * t676 - t936) * t821 * t869 + (((t658 * t986 + t661) * t793 + (-t736 - t827) * t926 - 0.2e1 * t881 - t661 + t676 * t985) * t861 + (t643 * t807 + (t963 + 0.1e1) * t700 * t990) * t866 + (t676 * t904 + t954 * t925) * t869) * qJ(3,1);
t874 = qJ(3,1) * t793 * t820;
t883 = -0.2e1 * t948;
t908 = t807 * t813;
t634 = (-t676 * t874 + (t676 * t990 - t661) * t908 + (-t792 * t914 + (qJ(3,1) * t883 - t827 + t891) * t813 - t746 * t793 * t882) * t700) * t861 + (-t700 * t874 + (t647 + t860) * t901 + t647 * t746) * t866 + (t676 * t746 - t936) * t869;
t640 = (0.2e1 * t807 * t664 - 0.2e1 * t735 * t676 + t925) * t700 * t725;
t697 = t700 ^ 2;
t758 = -mrSges(2,2) + t762;
t849 = -qJ(3,1) * mrSges(3,2) + t945;
t703 = -t813 * (t758 * pkin(5) - t849) + t726 * t807;
t944 = mrSges(3,3) * qJ(3,1);
t784 = 0.2e1 * t944;
t785 = -0.2e1 * t944;
t939 = t664 * t762;
t951 = pkin(1) * t758;
t895 = t703 * t640 + (-m(3) * (t827 + t830) + t785 + t877) * t634 + t981 * t631 + 0.2e1 * t676 * t939 + (t989 * t963 - ((t784 + t846 + t984) * t807 + t951) * t813 + t807 * t947 + t989) * t697;
t975 = t813 * t895;
t959 = m(3) * pkin(1);
t958 = m(3) * pkin(5);
t769 = mrSges(3,2) + t958;
t919 = t769 * t803;
t918 = t769 * t805;
t917 = t769 * t807;
t668 = t987 * t674;
t671 = t674 ^ 2;
t707 = (-mrSges(2,1) / 0.2e1 - mrSges(3,1) / 0.2e1) * pkin(5) + Ifges(3,4) / 0.2e1 + Ifges(2,5) / 0.2e1 + (-t958 / 0.2e1 - mrSges(3,2) / 0.2e1) * pkin(2);
t740 = 0.2e1 * t947;
t786 = 0.2e1 * t957;
t816 = mrSges(2,2) * pkin(5);
t838 = t887 * m(2) + (mrSges(3,2) + mrSges(2,3)) * t967 - Ifges(1,3) - t946;
t844 = -t856 + t982;
t900 = -(-(t781 + t786 - t844) * t789 - t740 * t809 + t756 * t885 - (t823 - t887) * m(3) + t781 - 0.2e1 * t987 * t912 + t838) * t638 - t701 * t632 + t629 * t919 - 0.4e1 * (t668 - t941 / 0.2e1) * t980 - 0.2e1 * (((t780 + t787 + t844) * t674 + t662 * t981) * t932 + t674 * (t662 * t769 + t707 * t674 + t698 * t953)) * t809 - ((-t760 * pkin(5) + t816 + t847) * t671 + (m(3) * t662 - t674 * t759) * t698 * t886) * t803 + 0.2e1 * (t668 - t941) * t698;
t669 = t988 * t675;
t672 = t675 ^ 2;
t843 = -t856 + t983;
t899 = -(-(t783 + t786 - t843) * t791 - t740 * t811 + t757 * t884 - (t825 - t887) * m(3) + t783 - 0.2e1 * t988 * t910 + t838) * t639 - t702 * t633 + t630 * t918 - 0.4e1 * (t669 - t940 / 0.2e1) * t979 - 0.2e1 * (((t782 + t787 + t843) * t675 + t663 * t981) * t929 + t675 * (t663 * t769 + t707 * t675 + t699 * t952)) * t811 - ((-t761 * pkin(5) + t816 + t848) * t672 + (m(3) * t663 - t675 * t759) * t699 * t886) * t805 + 0.2e1 * (t669 - t940) * t699;
t670 = t989 * t676;
t673 = t676 ^ 2;
t842 = -t856 + t984;
t898 = -(-(t785 + t786 - t842) * t793 - t740 * t813 + t758 * t883 - (t827 - t887) * m(3) + t785 - 0.2e1 * t989 * t908 + t838) * t640 - t703 * t634 + t631 * t917 - 0.4e1 * (t670 - t939 / 0.2e1) * t978 - 0.2e1 * (((t784 + t787 + t842) * t676 + t664 * t981) * t926 + t676 * (t664 * t769 + t707 * t676 + t700 * t951)) * t813 - ((-t762 * pkin(5) + t816 + t849) * t673 + (m(3) * t664 - t676 * t759) * t700 * t886) * t807 + 0.2e1 * (t670 - t939) * t700;
t894 = -m(3) * t629 + t981 * t632 - t638 * t919 - t671 * t760 + (t956 * t760 + (-t809 * t981 - t959) * t803) * t695;
t893 = -m(3) * t630 + t981 * t633 - t639 * t918 - t672 * t761 + (t955 * t761 + (-t811 * t981 - t959) * t805) * t696;
t892 = -m(3) * t631 + t981 * t634 - t640 * t917 - t673 * t762 + (t954 * t762 + (-t813 * t981 - t959) * t807) * t697;
t1 = [(t898 * t713 + (t892 * t687 - t693 * t975) * t828) * t725 + (t899 * t711 + (t893 * t685 - t691 * t976) * t826) * t724 + (t900 * t709 + (t894 * t683 - t689 * t977) * t824) * t723; (t898 * t712 + (t892 * t688 + t694 * t975) * t828) * t725 + (t899 * t710 + (t893 * t686 + t692 * t976) * t826) * t724 + (t900 * t708 + (t894 * t684 + t690 * t977) * t824) * t723; (t892 * t735 + t895 * t807) * t828 + (t893 * t734 + t896 * t805) * t826 + (t894 * t733 + t897 * t803) * t824;];
taucX  = t1;
