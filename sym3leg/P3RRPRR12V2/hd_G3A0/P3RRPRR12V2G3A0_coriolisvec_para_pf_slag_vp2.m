% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR12V2G3A0
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
% Datum: 2020-08-06 19:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:27:25
% EndTime: 2020-08-06 19:27:34
% DurationCPUTime: 8.71s
% Computational Cost: add. (79992->459), mult. (94779->718), div. (11403->6), fcn. (72603->18), ass. (0->292)
t818 = sin(qJ(2,3));
t825 = cos(qJ(1,3));
t819 = sin(qJ(1,3));
t835 = pkin(5) - pkin(6);
t897 = pkin(1) * t825 + t819 * t835;
t732 = qJ(3,3) * t825 + t897 * t818;
t775 = t818 * qJ(3,3);
t877 = t825 * t775;
t735 = 0.2e1 * t877 + t897;
t836 = pkin(2) + pkin(3);
t928 = (qJ(3,3) + t836) * (-qJ(3,3) + t836);
t871 = t818 * t928;
t982 = pkin(1) * qJ(3,3);
t742 = -t871 + t982;
t978 = pkin(1) * t818;
t762 = qJ(3,3) + t978;
t815 = legFrame(3,2);
t778 = sin(t815);
t781 = cos(t815);
t824 = cos(qJ(2,3));
t806 = t824 ^ 2;
t870 = t825 * t928;
t991 = -0.2e1 * t836;
t886 = qJ(3,3) * t991;
t934 = t781 * t836;
t937 = t778 * t836;
t966 = qJ(3,3) * t778;
t699 = (-t778 * t870 + t781 * t886) * t806 + (-t735 * t937 - t742 * t781) * t824 - t732 * t966 + t762 * t934;
t965 = qJ(3,3) * t781;
t700 = (t778 * t886 + t781 * t870) * t806 + (t735 * t934 - t742 * t778) * t824 + t732 * t965 + t762 * t937;
t872 = t806 * t928;
t914 = t835 * t825;
t919 = t824 * t836;
t723 = -t819 * t872 - ((0.2e1 * t775 + pkin(1)) * t819 - t914) * t919 - qJ(3,3) * (t762 * t819 - t818 * t914);
t832 = xDP(3);
t833 = xDP(2);
t834 = xDP(1);
t759 = t775 + pkin(1);
t998 = t759 + t919;
t745 = 0.1e1 / t998;
t839 = 0.1e1 / qJ(3,3);
t940 = t745 * t839;
t678 = (t699 * t833 + t700 * t834 + t723 * t832) * t940;
t736 = t877 + t897;
t918 = t825 * t836;
t924 = t818 * t836;
t708 = (-t778 * t918 - t965) * t806 + (-t736 * t778 + t781 * t924) * t824 + t781 * t762;
t709 = (t781 * t918 - t966) * t806 + (t736 * t781 + t778 * t924) * t824 + t778 * t762;
t943 = (t998 * t819 - t914) * t824;
t690 = (t708 * t833 + t709 * t834 - t832 * t943) * t940;
t958 = t690 * t836;
t673 = t678 - t958;
t820 = sin(qJ(2,2));
t827 = cos(qJ(1,2));
t821 = sin(qJ(1,2));
t896 = pkin(1) * t827 + t821 * t835;
t733 = qJ(3,2) * t827 + t896 * t820;
t776 = t820 * qJ(3,2);
t878 = t827 * t776;
t737 = 0.2e1 * t878 + t896;
t927 = (qJ(3,2) + t836) * (-qJ(3,2) + t836);
t868 = t820 * t927;
t983 = pkin(1) * qJ(3,2);
t743 = -t868 + t983;
t977 = pkin(1) * t820;
t763 = qJ(3,2) + t977;
t816 = legFrame(2,2);
t779 = sin(t816);
t782 = cos(t816);
t826 = cos(qJ(2,2));
t807 = t826 ^ 2;
t867 = t827 * t927;
t887 = qJ(3,2) * t991;
t933 = t782 * t836;
t936 = t779 * t836;
t968 = qJ(3,2) * t779;
t701 = (-t779 * t867 + t782 * t887) * t807 + (-t737 * t936 - t743 * t782) * t826 - t733 * t968 + t763 * t933;
t967 = qJ(3,2) * t782;
t702 = (t779 * t887 + t782 * t867) * t807 + (t737 * t933 - t743 * t779) * t826 + t733 * t967 + t763 * t936;
t869 = t807 * t927;
t910 = t836 * t826;
t913 = t835 * t827;
t724 = -t821 * t869 - ((0.2e1 * t776 + pkin(1)) * t821 - t913) * t910 - qJ(3,2) * (t763 * t821 - t820 * t913);
t760 = t776 + pkin(1);
t997 = t760 + t910;
t746 = 0.1e1 / t997;
t841 = 0.1e1 / qJ(3,2);
t939 = t746 * t841;
t679 = (t701 * t833 + t702 * t834 + t724 * t832) * t939;
t738 = t878 + t896;
t917 = t827 * t836;
t922 = t820 * t836;
t710 = (-t779 * t917 - t967) * t807 + (-t738 * t779 + t782 * t922) * t826 + t782 * t763;
t711 = (t782 * t917 - t968) * t807 + (t738 * t782 + t779 * t922) * t826 + t779 * t763;
t942 = (t997 * t821 - t913) * t826;
t691 = (t710 * t833 + t711 * t834 - t832 * t942) * t939;
t956 = t691 * t836;
t674 = t679 - t956;
t822 = sin(qJ(2,1));
t829 = cos(qJ(1,1));
t823 = sin(qJ(1,1));
t895 = pkin(1) * t829 + t823 * t835;
t734 = qJ(3,1) * t829 + t895 * t822;
t777 = t822 * qJ(3,1);
t879 = t829 * t777;
t739 = 0.2e1 * t879 + t895;
t926 = (qJ(3,1) + t836) * (-qJ(3,1) + t836);
t865 = t822 * t926;
t984 = pkin(1) * qJ(3,1);
t744 = -t865 + t984;
t976 = pkin(1) * t822;
t764 = qJ(3,1) + t976;
t817 = legFrame(1,2);
t780 = sin(t817);
t783 = cos(t817);
t828 = cos(qJ(2,1));
t808 = t828 ^ 2;
t864 = t829 * t926;
t888 = qJ(3,1) * t991;
t932 = t783 * t836;
t935 = t780 * t836;
t970 = qJ(3,1) * t780;
t703 = (-t780 * t864 + t783 * t888) * t808 + (-t739 * t935 - t744 * t783) * t828 - t734 * t970 + t764 * t932;
t969 = qJ(3,1) * t783;
t704 = (t780 * t888 + t783 * t864) * t808 + (t739 * t932 - t744 * t780) * t828 + t734 * t969 + t764 * t935;
t866 = t808 * t926;
t912 = t835 * t829;
t916 = t828 * t836;
t725 = -t823 * t866 - ((0.2e1 * t777 + pkin(1)) * t823 - t912) * t916 - qJ(3,1) * (t764 * t823 - t822 * t912);
t761 = t777 + pkin(1);
t996 = t761 + t916;
t747 = 0.1e1 / t996;
t843 = 0.1e1 / qJ(3,1);
t938 = t747 * t843;
t680 = (t703 * t833 + t704 * t834 + t725 * t832) * t938;
t740 = t879 + t895;
t915 = t829 * t836;
t920 = t822 * t836;
t712 = (-t780 * t915 - t969) * t808 + (-t740 * t780 + t783 * t920) * t828 + t783 * t764;
t713 = (t783 * t915 - t970) * t808 + (t740 * t783 + t780 * t920) * t828 + t780 * t764;
t941 = (t996 * t823 - t912) * t828;
t692 = (t712 * t833 + t713 * t834 - t832 * t941) * t938;
t954 = t692 * t836;
t672 = t680 - t954;
t1002 = m(3) * pkin(2) + mrSges(3,1);
t882 = pkin(2) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t1009 = qJ(3,1) * t1002 + t882;
t1008 = qJ(3,2) * t1002 + t882;
t1007 = qJ(3,3) * t1002 + t882;
t890 = 0.2e1 * pkin(1);
t1006 = 0.2e1 * t836;
t842 = qJ(3,1) ^ 2;
t846 = pkin(2) ^ 2;
t1005 = (t842 - t846) * m(3);
t840 = qJ(3,2) ^ 2;
t1004 = (t840 - t846) * m(3);
t838 = qJ(3,3) ^ 2;
t1003 = (t838 - t846) * m(3);
t774 = qJ(3,1) * m(3) + mrSges(3,3);
t772 = m(3) * qJ(3,3) + mrSges(3,3);
t773 = m(3) * qJ(3,2) + mrSges(3,3);
t957 = t690 * t838;
t663 = t673 * t836 - t957;
t670 = t835 * t673;
t720 = (-t825 * t832 + (t778 * t833 - t781 * t834) * t819) * t745;
t891 = -pkin(1) ^ 2 - pkin(5) ^ 2;
t995 = -0.2e1 * pkin(5);
t859 = -t891 + (t995 + pkin(6)) * pkin(6);
t755 = t838 + t859;
t810 = t836 ^ 2;
t883 = t720 * t982;
t900 = 0.2e1 * t883 + t670;
t911 = t835 * t836;
t925 = t818 * t835;
t950 = t720 * t835;
t951 = t720 * t818;
t952 = t720 * t806;
t875 = t720 * t925;
t961 = (t875 - t958) * t824;
t987 = t806 - 0.1e1;
t994 = -0.2e1 * t806;
t651 = ((-(t810 - 0.3e1 * t838) * t919 * t952 + (t835 * (-t1006 * t690 + t678) * qJ(3,3) + (-0.3e1 * (-t838 / 0.3e1 + t810) * t775 + (t838 - t810) * t890) * t720) * t806 + (-t925 * t957 + ((-0.4e1 * t883 - t670) * t818 - t720 * (0.3e1 * t838 + t859)) * t836) * t824 - qJ(3,3) * (t755 * t951 + t900)) * t720 + ((t663 * t836 + t871 * t950) * t824 + t663 * pkin(1) + (t663 * t818 + (t994 + 0.1e1) * t720 * t911) * qJ(3,3)) * t690 + ((pkin(1) * t690 - t961) * t836 + (t690 * t924 + t987 * t950) * qJ(3,3)) * t678) * t940;
t863 = t824 * qJ(3,3) * t690;
t654 = (-(t835 * t863 + t900 * t818 + (t1006 * t759 * t824 + t755 + t872) * t720) * t824 * t720 + (-qJ(3,3) * t806 * t950 + (t673 + t875) * t919 + t673 * t759) * t690 + (t690 * t759 - t961) * t678) * t940;
t660 = (0.2e1 * t673 * t818 + 0.2e1 * t863 + t950) * t720 * t745;
t684 = t1007 * t690;
t687 = t690 ^ 2;
t771 = mrSges(2,1) + t1002;
t748 = pkin(2) * mrSges(3,2) + pkin(5) * t771 - Ifges(3,4) - Ifges(2,5);
t768 = -mrSges(2,2) + t772;
t974 = Ifges(3,6) - Ifges(2,6);
t856 = -qJ(3,3) * mrSges(3,2) + t974;
t726 = -t824 * (pkin(5) * t768 - t856) + t818 * t748;
t989 = m(3) * pkin(5);
t741 = (-mrSges(2,1) / 0.2e1 - mrSges(3,1) / 0.2e1) * pkin(5) + Ifges(3,4) / 0.2e1 + Ifges(2,5) / 0.2e1 + (-t989 / 0.2e1 - mrSges(3,2) / 0.2e1) * pkin(2);
t758 = t771 * t890;
t784 = mrSges(3,2) + t989;
t971 = mrSges(3,3) * qJ(3,3);
t798 = 0.2e1 * t971;
t799 = -0.2e1 * t971;
t988 = mrSges(3,1) * pkin(2);
t804 = 0.2e1 * t988;
t805 = -0.2e1 * t988;
t831 = mrSges(2,2) * pkin(5);
t975 = Ifges(2,1) + Ifges(3,1);
t851 = t891 * m(2) + (mrSges(3,2) + mrSges(2,3)) * t995 - Ifges(1,3) - t975;
t861 = Ifges(2,2) + Ifges(3,3) - t975;
t854 = -t861 + t1003;
t931 = t784 * t818;
t964 = t678 * t772;
t981 = pkin(1) * t768;
t909 = (-(t799 + t804 - t854) * t806 - 0.2e1 * t768 * t978 - (t838 - t891) * m(3) + t799 + t851 + (-0.2e1 * t1007 * t818 - t758) * t824) * t660 + t726 * t654 - t651 * t931 + 0.4e1 * (t684 - t964 / 0.2e1) * t952 + 0.2e1 * (((t798 + t805 + t854) * t690 + t678 * t1002) * t951 + t690 * (t678 * t784 + t690 * t741 + t720 * t981)) * t824 + ((-t772 * pkin(5) + t831 + t856) * t687 + (m(3) * t678 - t690 * t771) * t720 * t890) * t818 - 0.2e1 * t720 * (t684 - t964);
t1001 = t819 * t909;
t955 = t691 * t840;
t664 = t674 * t836 - t955;
t671 = t835 * t674;
t721 = (-t827 * t832 + (t779 * t833 - t782 * t834) * t821) * t746;
t756 = t840 + t859;
t884 = t721 * t983;
t899 = 0.2e1 * t884 + t671;
t923 = t820 * t835;
t947 = t721 * t835;
t948 = t721 * t820;
t949 = t721 * t807;
t874 = t721 * t923;
t960 = (t874 - t956) * t826;
t986 = t807 - 0.1e1;
t993 = -0.2e1 * t807;
t652 = ((-(t810 - 0.3e1 * t840) * t910 * t949 + (t835 * (-t1006 * t691 + t679) * qJ(3,2) + (-0.3e1 * (-t840 / 0.3e1 + t810) * t776 + (t840 - t810) * t890) * t721) * t807 + (-t923 * t955 + ((-0.4e1 * t884 - t671) * t820 - t721 * (0.3e1 * t840 + t859)) * t836) * t826 - qJ(3,2) * (t756 * t948 + t899)) * t721 + ((t664 * t836 + t868 * t947) * t826 + t664 * pkin(1) + (t664 * t820 + (t993 + 0.1e1) * t721 * t911) * qJ(3,2)) * t691 + ((pkin(1) * t691 - t960) * t836 + (t691 * t922 + t986 * t947) * qJ(3,2)) * t679) * t939;
t862 = t826 * qJ(3,2) * t691;
t655 = (-(t835 * t862 + t899 * t820 + (t1006 * t760 * t826 + t756 + t869) * t721) * t826 * t721 + (-qJ(3,2) * t807 * t947 + (t674 + t874) * t910 + t674 * t760) * t691 + (t691 * t760 - t960) * t679) * t939;
t661 = (0.2e1 * t674 * t820 + 0.2e1 * t862 + t947) * t721 * t746;
t685 = t1008 * t691;
t688 = t691 ^ 2;
t769 = -mrSges(2,2) + t773;
t857 = -qJ(3,2) * mrSges(3,2) + t974;
t727 = -t826 * (pkin(5) * t769 - t857) + t820 * t748;
t972 = mrSges(3,3) * qJ(3,2);
t800 = 0.2e1 * t972;
t801 = -0.2e1 * t972;
t853 = -t861 + t1004;
t930 = t784 * t820;
t963 = t679 * t773;
t980 = pkin(1) * t769;
t908 = (-(t801 + t804 - t853) * t807 - 0.2e1 * t769 * t977 - (t840 - t891) * m(3) + t801 + t851 + (-0.2e1 * t1008 * t820 - t758) * t826) * t661 + t727 * t655 - t652 * t930 + 0.4e1 * (t685 - t963 / 0.2e1) * t949 + 0.2e1 * (((t800 + t805 + t853) * t691 + t679 * t1002) * t948 + t691 * (t679 * t784 + t691 * t741 + t721 * t980)) * t826 + ((-t773 * pkin(5) + t831 + t857) * t688 + (m(3) * t679 - t691 * t771) * t721 * t890) * t820 - 0.2e1 * t721 * (t685 - t963);
t1000 = t821 * t908;
t953 = t692 * t842;
t665 = t672 * t836 - t953;
t669 = t835 * t672;
t722 = (-t829 * t832 + (t780 * t833 - t783 * t834) * t823) * t747;
t757 = t842 + t859;
t885 = t722 * t984;
t898 = 0.2e1 * t885 + t669;
t921 = t822 * t835;
t944 = t722 * t835;
t945 = t722 * t822;
t946 = t722 * t808;
t873 = t722 * t921;
t959 = (t873 - t954) * t828;
t985 = t808 - 0.1e1;
t992 = -0.2e1 * t808;
t653 = ((-(t810 - 0.3e1 * t842) * t916 * t946 + (t835 * (-t1006 * t692 + t680) * qJ(3,1) + (-0.3e1 * (-t842 / 0.3e1 + t810) * t777 + (t842 - t810) * t890) * t722) * t808 + (-t921 * t953 + ((-0.4e1 * t885 - t669) * t822 - t722 * (0.3e1 * t842 + t859)) * t836) * t828 - qJ(3,1) * (t757 * t945 + t898)) * t722 + ((t665 * t836 + t865 * t944) * t828 + t665 * pkin(1) + (t665 * t822 + (t992 + 0.1e1) * t722 * t911) * qJ(3,1)) * t692 + ((pkin(1) * t692 - t959) * t836 + (t692 * t920 + t985 * t944) * qJ(3,1)) * t680) * t938;
t876 = t692 * t828 * qJ(3,1);
t656 = (-(t835 * t876 + t898 * t822 + (t1006 * t761 * t828 + t757 + t866) * t722) * t828 * t722 + (-qJ(3,1) * t808 * t944 + (t672 + t873) * t916 + t672 * t761) * t692 + (t692 * t761 - t959) * t680) * t938;
t662 = (0.2e1 * t672 * t822 + 0.2e1 * t876 + t944) * t722 * t747;
t686 = t1009 * t692;
t689 = t692 ^ 2;
t770 = -mrSges(2,2) + t774;
t858 = -qJ(3,1) * mrSges(3,2) + t974;
t728 = -t828 * (pkin(5) * t770 - t858) + t822 * t748;
t973 = mrSges(3,3) * qJ(3,1);
t802 = 0.2e1 * t973;
t803 = -0.2e1 * t973;
t852 = -t861 + t1005;
t929 = t784 * t822;
t962 = t680 * t774;
t979 = pkin(1) * t770;
t907 = (-(t803 + t804 - t852) * t808 - 0.2e1 * t770 * t976 - (t842 - t891) * m(3) + t803 + t851 + (-0.2e1 * t1009 * t822 - t758) * t828) * t662 + t728 * t656 - t653 * t929 + 0.4e1 * (t686 - t962 / 0.2e1) * t946 + 0.2e1 * (((t802 + t805 + t852) * t692 + t680 * t1002) * t945 + t692 * (t680 * t784 + t692 * t741 + t722 * t979)) * t828 + ((-t774 * pkin(5) + t831 + t858) * t689 + (m(3) * t680 - t692 * t771) * t722 * t890) * t822 - 0.2e1 * t722 * (t686 - t962);
t999 = t823 * t907;
t990 = m(3) * pkin(1);
t717 = t720 ^ 2;
t855 = t805 - t861;
t881 = t805 - Ifges(3,2) - Ifges(2,3);
t906 = t726 * t660 + (-m(3) * (t838 + t846) + t799 + t881) * t654 + t1002 * t651 + 0.2e1 * t690 * t964 + (t1007 * t994 - ((t798 + t855 + t1003) * t818 + t981) * t824 + t771 * t978 + t1007) * t717;
t718 = t721 ^ 2;
t905 = t727 * t661 + (-m(3) * (t840 + t846) + t801 + t881) * t655 + t1002 * t652 + 0.2e1 * t691 * t963 + (t1008 * t993 - ((t800 + t855 + t1004) * t820 + t980) * t826 + t771 * t977 + t1008) * t718;
t719 = t722 ^ 2;
t904 = t728 * t662 + (-m(3) * (t842 + t846) + t803 + t881) * t656 + t1002 * t653 + 0.2e1 * t692 * t962 + (t1009 * t992 - ((t802 + t855 + t1005) * t822 + t979) * t828 + t771 * t976 + t1009) * t719;
t903 = -m(3) * t651 + t654 * t1002 - t660 * t931 - t687 * t772 + (t987 * t772 + (-t1002 * t824 - t990) * t818) * t717;
t902 = -m(3) * t652 + t655 * t1002 - t661 * t930 - t688 * t773 + (t986 * t773 + (-t1002 * t826 - t990) * t820) * t718;
t901 = -m(3) * t653 + t656 * t1002 - t662 * t929 - t689 * t774 + (t985 * t774 + (-t1002 * t828 - t990) * t822) * t719;
t1 = [(-t783 * t999 + (t901 * t704 + t904 * t713) * t843) * t747 + (-t782 * t1000 + (t902 * t702 + t905 * t711) * t841) * t746 + (-t781 * t1001 + (t903 * t700 + t906 * t709) * t839) * t745; (t780 * t999 + (t901 * t703 + t904 * t712) * t843) * t747 + (t779 * t1000 + (t902 * t701 + t905 * t710) * t841) * t746 + (t778 * t1001 + (t903 * t699 + t906 * t708) * t839) * t745; (-t907 * t829 + (t901 * t725 - t904 * t941) * t843) * t747 + (-t908 * t827 + (t902 * t724 - t905 * t942) * t841) * t746 + (-t909 * t825 + (t903 * t723 - t906 * t943) * t839) * t745;];
taucX  = t1;
