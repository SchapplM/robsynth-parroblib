% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR12V1G1A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
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
% Datum: 2020-08-06 19:02
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:01:47
% EndTime: 2020-08-06 19:01:51
% DurationCPUTime: 4.02s
% Computational Cost: add. (19812->269), mult. (32112->437), div. (4653->6), fcn. (29763->18), ass. (0->215)
t695 = sin(qJ(2,3));
t701 = cos(qJ(2,3));
t711 = pkin(1) + pkin(2);
t761 = t711 * t701;
t770 = t695 * qJ(3,3);
t652 = t761 + t770;
t702 = cos(qJ(1,3));
t672 = t702 * pkin(4);
t696 = sin(qJ(1,3));
t616 = t652 * t696 + t672;
t790 = t696 * pkin(4);
t619 = t652 * t702 - t790;
t689 = legFrame(3,3);
t662 = sin(t689);
t665 = cos(t689);
t610 = -t662 * t616 + t665 * t619;
t613 = t616 * t665 + t662 * t619;
t708 = xDP(3);
t709 = xDP(2);
t710 = xDP(1);
t713 = 0.1e1 / qJ(3,3);
t764 = t711 * t695;
t725 = t701 * qJ(3,3) - t764;
t586 = (t610 * t710 + t613 * t709 - t708 * t725) * t713;
t625 = -t662 * t696 + t665 * t702;
t637 = t696 * t770 + t672;
t638 = t702 * t770 - t790;
t598 = t625 * t761 - t662 * t637 + t638 * t665;
t626 = t662 * t702 + t665 * t696;
t599 = t626 * t761 + t637 * t665 + t662 * t638;
t643 = 0.1e1 / t652;
t583 = (t695 * t708 + (t598 * t710 + t599 * t709) * t701 * t643) * t713;
t812 = t711 * t583;
t818 = t586 - t812;
t821 = t695 * t818;
t697 = sin(qJ(2,2));
t703 = cos(qJ(2,2));
t760 = t711 * t703;
t768 = t697 * qJ(3,2);
t653 = t760 + t768;
t704 = cos(qJ(1,2));
t673 = t704 * pkin(4);
t698 = sin(qJ(1,2));
t617 = t653 * t698 + t673;
t789 = t698 * pkin(4);
t620 = t653 * t704 - t789;
t690 = legFrame(2,3);
t663 = sin(t690);
t666 = cos(t690);
t611 = -t663 * t617 + t666 * t620;
t614 = t617 * t666 + t663 * t620;
t715 = 0.1e1 / qJ(3,2);
t763 = t711 * t697;
t726 = t703 * qJ(3,2) - t763;
t587 = (t611 * t710 + t614 * t709 - t708 * t726) * t715;
t627 = -t663 * t698 + t666 * t704;
t639 = t698 * t768 + t673;
t640 = t704 * t768 - t789;
t600 = t627 * t760 - t663 * t639 + t640 * t666;
t628 = t663 * t704 + t666 * t698;
t601 = t628 * t760 + t639 * t666 + t663 * t640;
t644 = 0.1e1 / t653;
t584 = (t697 * t708 + (t600 * t710 + t601 * t709) * t703 * t644) * t715;
t811 = t711 * t584;
t817 = t587 - t811;
t820 = t697 * t817;
t699 = sin(qJ(2,1));
t705 = cos(qJ(2,1));
t759 = t711 * t705;
t766 = t699 * qJ(3,1);
t654 = t759 + t766;
t706 = cos(qJ(1,1));
t674 = t706 * pkin(4);
t700 = sin(qJ(1,1));
t618 = t654 * t700 + t674;
t788 = t700 * pkin(4);
t621 = t654 * t706 - t788;
t691 = legFrame(1,3);
t664 = sin(t691);
t667 = cos(t691);
t612 = -t664 * t618 + t667 * t621;
t615 = t618 * t667 + t664 * t621;
t717 = 0.1e1 / qJ(3,1);
t762 = t711 * t699;
t727 = t705 * qJ(3,1) - t762;
t588 = (t612 * t710 + t615 * t709 - t708 * t727) * t717;
t629 = -t664 * t700 + t667 * t706;
t641 = t700 * t766 + t674;
t642 = t706 * t766 - t788;
t602 = t629 * t759 - t664 * t641 + t642 * t667;
t630 = t664 * t706 + t667 * t700;
t603 = t630 * t759 + t641 * t667 + t664 * t642;
t645 = 0.1e1 / t654;
t585 = (t699 * t708 + (t602 * t710 + t603 * t709) * t705 * t645) * t717;
t810 = t711 * t585;
t816 = t588 - t810;
t819 = t699 * t816;
t815 = 0.2e1 * t583;
t814 = 0.2e1 * t584;
t813 = 0.2e1 * t585;
t712 = qJ(3,3) ^ 2;
t801 = -0.2e1 * qJ(3,3);
t747 = mrSges(3,3) * t801;
t809 = -m(3) * t712 + t747;
t714 = qJ(3,2) ^ 2;
t802 = -0.2e1 * qJ(3,2);
t748 = mrSges(3,3) * t802;
t808 = -m(3) * t714 + t748;
t716 = qJ(3,1) ^ 2;
t803 = -0.2e1 * qJ(3,1);
t749 = mrSges(3,3) * t803;
t807 = -m(3) * t716 + t749;
t595 = (t625 * t709 - t626 * t710) * t643;
t685 = t701 ^ 2;
t718 = pkin(4) ^ 2;
t806 = pkin(4) * t821 + (-t712 - t718 - (qJ(3,3) + t711) * (-qJ(3,3) + t711) * t685) * t595;
t596 = (t627 * t709 - t628 * t710) * t644;
t686 = t703 ^ 2;
t805 = pkin(4) * t820 + (-t714 - t718 - (qJ(3,2) + t711) * (-qJ(3,2) + t711) * t686) * t596;
t597 = (t629 * t709 - t630 * t710) * t645;
t687 = t705 ^ 2;
t804 = pkin(4) * t819 + (-t716 - t718 - (qJ(3,1) + t711) * (-qJ(3,1) + t711) * t687) * t597;
t800 = 0.2e1 * t586;
t799 = 0.2e1 * t587;
t798 = 0.2e1 * t588;
t668 = m(3) * pkin(1) + mrSges(3,1);
t742 = pkin(1) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t646 = t668 * qJ(3,3) + t742;
t797 = -0.2e1 * t646;
t647 = t668 * qJ(3,2) + t742;
t796 = -0.2e1 * t647;
t648 = t668 * qJ(3,1) + t742;
t795 = -0.2e1 * t648;
t794 = mrSges(3,1) * pkin(1);
t793 = pkin(4) * t583;
t792 = pkin(4) * t584;
t791 = pkin(4) * t585;
t787 = -Ifges(2,1) - Ifges(3,1);
t786 = Ifges(2,6) - Ifges(3,6);
t785 = mrSges(3,2) * t695;
t784 = mrSges(3,2) * t697;
t783 = mrSges(3,2) * t699;
t659 = m(3) * qJ(3,3) + mrSges(3,3);
t782 = t586 * t659;
t660 = m(3) * qJ(3,2) + mrSges(3,3);
t781 = t587 * t660;
t661 = m(3) * qJ(3,1) + mrSges(3,3);
t780 = t588 * t661;
t779 = t595 * t685;
t778 = t595 * t695;
t777 = t596 * t686;
t776 = t596 * t697;
t775 = t597 * t687;
t774 = t597 * t699;
t769 = t695 * t701;
t767 = t697 * t703;
t765 = t699 * t705;
t589 = pkin(4) * t778;
t568 = t589 + t812;
t744 = (t595 * t764 - t793 / 0.2e1) * t801;
t544 = ((t685 * t744 + t806 * t701) * t595 + (-(t589 - t818) * t761 + (pkin(4) * t779 + t821) * qJ(3,3)) * t583 + (t568 * t701 + t583 * t770) * t586) * t643 * t713;
t720 = pkin(1) ^ 2;
t737 = -t720 + (-0.2e1 * pkin(1) - pkin(2)) * pkin(2);
t547 = (((-t712 + t737) * t583 + t711 * t586) * t583 + t568 * t586 + (t701 * t744 + t725 * t793 + t806) * t595) * t713;
t553 = (-pkin(4) * t595 + t695 * t800 + t725 * t815) * t595 * t643;
t571 = t646 * t583;
t580 = t583 ^ 2;
t655 = mrSges(3,2) * qJ(3,3) + t786;
t658 = pkin(1) * mrSges(3,2) - Ifges(3,4) - Ifges(2,5);
t622 = -t655 * t701 + t695 * t658;
t724 = m(3) * t720 + Ifges(2,2) + Ifges(3,3) + t787 + 0.2e1 * t794;
t631 = t724 + t809;
t743 = -Ifges(1,3) + t787;
t758 = (-t631 * t685 + t769 * t797 + t743 + t809) * t553 + t622 * t544 - t547 * t785 + 0.4e1 * (t571 - t782 / 0.2e1) * t779 + (0.2e1 * (-t631 * t583 + t586 * t668) * t778 + (mrSges(3,2) * t800 - t658 * t583) * t583) * t701 - t580 * t655 * t695 - 0.2e1 * (t571 - t782) * t595;
t590 = pkin(4) * t776;
t569 = t590 + t811;
t745 = (t596 * t763 - t792 / 0.2e1) * t802;
t545 = ((t686 * t745 + t805 * t703) * t596 + (-(t590 - t817) * t760 + (pkin(4) * t777 + t820) * qJ(3,2)) * t584 + (t569 * t703 + t584 * t768) * t587) * t644 * t715;
t548 = (((-t714 + t737) * t584 + t711 * t587) * t584 + t569 * t587 + (t703 * t745 + t726 * t792 + t805) * t596) * t715;
t554 = (-pkin(4) * t596 + t697 * t799 + t726 * t814) * t596 * t644;
t572 = t647 * t584;
t581 = t584 ^ 2;
t656 = mrSges(3,2) * qJ(3,2) + t786;
t623 = -t656 * t703 + t697 * t658;
t632 = t724 + t808;
t757 = (-t632 * t686 + t767 * t796 + t743 + t808) * t554 + t623 * t545 - t548 * t784 + 0.4e1 * (t572 - t781 / 0.2e1) * t777 + (0.2e1 * (-t632 * t584 + t587 * t668) * t776 + (mrSges(3,2) * t799 - t658 * t584) * t584) * t703 - t581 * t656 * t697 - 0.2e1 * (t572 - t781) * t596;
t591 = pkin(4) * t774;
t570 = t591 + t810;
t746 = (t597 * t762 - t791 / 0.2e1) * t803;
t546 = ((t687 * t746 + t804 * t705) * t597 + (-(t591 - t816) * t759 + (pkin(4) * t775 + t819) * qJ(3,1)) * t585 + (t570 * t705 + t585 * t766) * t588) * t645 * t717;
t549 = (((-t716 + t737) * t585 + t711 * t588) * t585 + t570 * t588 + (t705 * t746 + t727 * t791 + t804) * t597) * t717;
t555 = (-pkin(4) * t597 + t699 * t798 + t727 * t813) * t597 * t645;
t573 = t648 * t585;
t582 = t585 ^ 2;
t657 = qJ(3,1) * mrSges(3,2) + t786;
t624 = -t657 * t705 + t699 * t658;
t633 = t724 + t807;
t756 = (-t633 * t687 + t765 * t795 + t743 + t807) * t555 + t624 * t546 - t549 * t783 + 0.4e1 * (t573 - t780 / 0.2e1) * t775 + (0.2e1 * (-t633 * t585 + t588 * t668) * t774 + (mrSges(3,2) * t798 - t658 * t585) * t585) * t705 - t582 * t657 * t699 - 0.2e1 * (t573 - t780) * t597;
t592 = t595 ^ 2;
t741 = -0.2e1 * t794 - Ifges(3,2) - Ifges(2,3);
t755 = t622 * t553 + (-(t712 + t720) * m(3) + t747 + t741) * t544 + t668 * t547 + t782 * t815 + (t631 * t769 + t685 * t797 + t646) * t592;
t593 = t596 ^ 2;
t754 = t623 * t554 + (-(t714 + t720) * m(3) + t748 + t741) * t545 + t668 * t548 + t781 * t814 + (t632 * t767 + t686 * t796 + t647) * t593;
t594 = t597 ^ 2;
t753 = t624 * t555 + (-(t716 + t720) * m(3) + t749 + t741) * t546 + t668 * t549 + t780 * t813 + (t633 * t765 + t687 * t795 + t648) * t594;
t752 = -m(3) * t547 + t668 * t544 - t553 * t785 - t580 * t659 + (t659 * t685 - t668 * t769 - t659) * t592;
t751 = -m(3) * t548 + t668 * t545 - t554 * t784 - t581 * t660 + (t660 * t686 - t668 * t767 - t660) * t593;
t750 = -m(3) * t549 + t668 * t546 - t555 * t783 - t582 * t661 + (t661 * t687 - t668 * t765 - t661) * t594;
t736 = t752 * t713;
t735 = t751 * t715;
t734 = t750 * t717;
t733 = t755 * t713 * t701;
t732 = t754 * t715 * t703;
t731 = t753 * t717 * t705;
t1 = [t612 * t734 + t611 * t735 + t610 * t736 + (t602 * t731 - t756 * t630) * t645 + (t600 * t732 - t757 * t628) * t644 + (t598 * t733 - t758 * t626) * t643; t615 * t734 + t614 * t735 + t613 * t736 + (t603 * t731 + t756 * t629) * t645 + (t601 * t732 + t757 * t627) * t644 + (t599 * t733 + t758 * t625) * t643; (t753 * t699 - t727 * t750) * t717 + (t754 * t697 - t726 * t751) * t715 + (t755 * t695 - t725 * t752) * t713;];
taucX  = t1;
