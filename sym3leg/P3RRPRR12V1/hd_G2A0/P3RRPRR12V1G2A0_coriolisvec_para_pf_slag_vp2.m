% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR12V1G2A0
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
% Datum: 2020-08-06 19:07
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:05:45
% EndTime: 2020-08-06 19:05:49
% DurationCPUTime: 3.52s
% Computational Cost: add. (27090->275), mult. (41871->447), div. (6471->6), fcn. (35190->18), ass. (0->212)
t667 = sin(qJ(1,3));
t666 = sin(qJ(2,3));
t755 = qJ(3,3) * t666;
t673 = cos(qJ(1,3));
t766 = pkin(4) * t673;
t623 = t667 * t755 + t766;
t663 = legFrame(3,2);
t645 = sin(t663);
t648 = cos(t663);
t672 = cos(qJ(2,3));
t656 = t672 ^ 2;
t682 = pkin(1) + pkin(2);
t716 = t682 * t666;
t725 = t667 * t682;
t732 = t645 * qJ(3,3);
t587 = (t648 * t725 - t732) * t656 + (t623 * t648 + t645 * t716) * t672 + t732;
t729 = t648 * qJ(3,3);
t590 = (-t645 * t725 - t729) * t656 + (-t623 * t645 + t648 * t716) * t672 + t729;
t679 = xDP(3);
t680 = xDP(2);
t681 = xDP(1);
t719 = t672 * t682;
t635 = t719 + t755;
t626 = 0.1e1 / t635;
t684 = 0.1e1 / qJ(3,3);
t735 = t626 * t684;
t698 = -pkin(4) * t667 + t635 * t673;
t738 = t698 * t672;
t572 = (t587 * t681 + t590 * t680 + t679 * t738) * t735;
t786 = t572 * t682;
t669 = sin(qJ(1,2));
t668 = sin(qJ(2,2));
t757 = qJ(3,2) * t668;
t675 = cos(qJ(1,2));
t765 = pkin(4) * t675;
t624 = t669 * t757 + t765;
t664 = legFrame(2,2);
t646 = sin(t664);
t649 = cos(t664);
t674 = cos(qJ(2,2));
t657 = t674 ^ 2;
t715 = t682 * t668;
t723 = t669 * t682;
t731 = t646 * qJ(3,2);
t588 = (t649 * t723 - t731) * t657 + (t624 * t649 + t646 * t715) * t674 + t731;
t728 = t649 * qJ(3,2);
t591 = (-t646 * t723 - t728) * t657 + (-t624 * t646 + t649 * t715) * t674 + t728;
t718 = t674 * t682;
t636 = t718 + t757;
t627 = 0.1e1 / t636;
t686 = 0.1e1 / qJ(3,2);
t734 = t627 * t686;
t697 = -pkin(4) * t669 + t636 * t675;
t737 = t697 * t674;
t573 = (t588 * t681 + t591 * t680 + t679 * t737) * t734;
t785 = t573 * t682;
t671 = sin(qJ(1,1));
t670 = sin(qJ(2,1));
t758 = qJ(3,1) * t670;
t677 = cos(qJ(1,1));
t764 = pkin(4) * t677;
t625 = t671 * t758 + t764;
t665 = legFrame(1,2);
t647 = sin(t665);
t650 = cos(t665);
t676 = cos(qJ(2,1));
t658 = t676 ^ 2;
t720 = t671 * t682;
t721 = t670 * t682;
t730 = t647 * qJ(3,1);
t589 = (t650 * t720 - t730) * t658 + (t625 * t650 + t647 * t721) * t676 + t730;
t727 = t650 * qJ(3,1);
t592 = (-t647 * t720 - t727) * t658 + (-t625 * t647 + t650 * t721) * t676 + t727;
t717 = t676 * t682;
t637 = t717 + t758;
t628 = 0.1e1 / t637;
t688 = 0.1e1 / qJ(3,1);
t733 = t628 * t688;
t696 = -pkin(4) * t671 + t637 * t677;
t736 = t696 * t676;
t574 = (t589 * t681 + t592 * t680 + t679 * t736) * t733;
t784 = t574 * t682;
t754 = qJ(3,3) * t672;
t693 = -t716 + t754;
t783 = t693 * t572;
t756 = qJ(3,2) * t674;
t694 = -t715 + t756;
t782 = t694 * t573;
t695 = qJ(3,1) * t676 - t721;
t781 = t695 * t574;
t584 = (-t667 * t679 + (-t645 * t680 + t648 * t681) * t673) * t626;
t683 = qJ(3,3) ^ 2;
t689 = pkin(4) ^ 2;
t750 = t572 * qJ(3,3);
t608 = t635 * t667 + t766;
t599 = t608 * t648 - t645 * t693;
t600 = -t608 * t645 - t648 * t693;
t575 = (t599 * t681 + t600 * t680 + t679 * t698) * t684;
t551 = -t575 + t786;
t753 = t551 * t666;
t545 = (-t672 * t750 + t753) * pkin(4) + ((qJ(3,3) + t682) * (-qJ(3,3) + t682) * t656 + 0.2e1 * t716 * t754 + t683 + t689) * t584;
t743 = t584 * t666;
t578 = pkin(4) * t743;
t557 = t578 + t786;
t744 = t584 * t656;
t533 = (-t545 * t672 * t584 + (-(t578 + t551) * t719 + (pkin(4) * t744 - t753) * qJ(3,3)) * t572 + (t557 * t672 + t666 * t750) * t575) * t735;
t691 = (pkin(1) ^ 2);
t705 = -t691 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t536 = (((-t683 + t705) * t572 + t682 * t575) * t572 + t557 * t575 + (pkin(4) * t783 - t545) * t584) * t684;
t773 = 0.2e1 * t575;
t542 = (-pkin(4) * t584 + t666 * t773 + 0.2e1 * t783) * t584 * t626;
t651 = m(3) * pkin(1) + mrSges(3,1);
t707 = pkin(1) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t629 = qJ(3,3) * t651 + t707;
t560 = t629 * t572;
t569 = t572 ^ 2;
t762 = Ifges(2,6) - Ifges(3,6);
t638 = mrSges(3,2) * qJ(3,3) + t762;
t641 = pkin(1) * mrSges(3,2) - Ifges(3,4) - Ifges(2,5);
t617 = -t638 * t672 + t641 * t666;
t763 = (-Ifges(2,1) - Ifges(3,1));
t767 = mrSges(3,1) * pkin(1);
t692 = m(3) * t691 + Ifges(2,2) + Ifges(3,3) + t763 + 2 * t767;
t774 = -2 * mrSges(3,3);
t709 = qJ(3,3) * t774;
t777 = -m(3) * t683 + t709;
t620 = t692 + t777;
t708 = -Ifges(1,3) + t763;
t726 = t666 * t672;
t642 = m(3) * qJ(3,3) + mrSges(3,3);
t747 = t575 * t642;
t761 = mrSges(3,2) * t666;
t770 = -0.2e1 * t629;
t714 = (-t620 * t656 + t726 * t770 + t708 + t777) * t542 + t617 * t533 - t536 * t761 + 0.4e1 * (t560 - t747 / 0.2e1) * t744 + (0.2e1 * (-t572 * t620 + t575 * t651) * t743 + (mrSges(3,2) * t773 - t572 * t641) * t572) * t672 - t569 * t638 * t666 - 0.2e1 * (t560 - t747) * t584;
t780 = t673 * t714;
t585 = (-t669 * t679 + (-t646 * t680 + t649 * t681) * t675) * t627;
t685 = qJ(3,2) ^ 2;
t749 = t573 * qJ(3,2);
t609 = t636 * t669 + t765;
t601 = t609 * t649 - t646 * t694;
t602 = -t609 * t646 - t649 * t694;
t576 = (t601 * t681 + t602 * t680 + t679 * t697) * t686;
t552 = -t576 + t785;
t752 = t552 * t668;
t546 = (-t674 * t749 + t752) * pkin(4) + ((qJ(3,2) + t682) * (-qJ(3,2) + t682) * t657 + 0.2e1 * t715 * t756 + t685 + t689) * t585;
t741 = t585 * t668;
t579 = pkin(4) * t741;
t558 = t579 + t785;
t742 = t585 * t657;
t534 = (-t546 * t674 * t585 + (-(t579 + t552) * t718 + (pkin(4) * t742 - t752) * qJ(3,2)) * t573 + (t558 * t674 + t668 * t749) * t576) * t734;
t537 = (((-t685 + t705) * t573 + t682 * t576) * t573 + t558 * t576 + (pkin(4) * t782 - t546) * t585) * t686;
t772 = 0.2e1 * t576;
t543 = (-pkin(4) * t585 + t668 * t772 + 0.2e1 * t782) * t585 * t627;
t630 = qJ(3,2) * t651 + t707;
t561 = t630 * t573;
t570 = t573 ^ 2;
t639 = mrSges(3,2) * qJ(3,2) + t762;
t618 = -t639 * t674 + t641 * t668;
t710 = qJ(3,2) * t774;
t776 = -m(3) * t685 + t710;
t621 = t692 + t776;
t724 = t668 * t674;
t643 = m(3) * qJ(3,2) + mrSges(3,3);
t746 = t576 * t643;
t760 = mrSges(3,2) * t668;
t769 = -0.2e1 * t630;
t713 = (-t621 * t657 + t724 * t769 + t708 + t776) * t543 + t618 * t534 - t537 * t760 + 0.4e1 * (t561 - t746 / 0.2e1) * t742 + (0.2e1 * (-t573 * t621 + t576 * t651) * t741 + (mrSges(3,2) * t772 - t573 * t641) * t573) * t674 - t570 * t639 * t668 - 0.2e1 * (t561 - t746) * t585;
t779 = t675 * t713;
t586 = (-t671 * t679 + (-t647 * t680 + t650 * t681) * t677) * t628;
t687 = qJ(3,1) ^ 2;
t748 = t574 * qJ(3,1);
t610 = t637 * t671 + t764;
t603 = t610 * t650 - t647 * t695;
t604 = -t610 * t647 - t650 * t695;
t577 = (t603 * t681 + t604 * t680 + t679 * t696) * t688;
t553 = -t577 + t784;
t751 = t553 * t670;
t547 = (-t676 * t748 + t751) * pkin(4) + ((qJ(3,1) + t682) * (-qJ(3,1) + t682) * t658 + 0.2e1 * t717 * t758 + t687 + t689) * t586;
t739 = t586 * t670;
t580 = pkin(4) * t739;
t559 = t580 + t784;
t740 = t586 * t658;
t535 = (-t547 * t676 * t586 + (-(t580 + t553) * t717 + (pkin(4) * t740 - t751) * qJ(3,1)) * t574 + (t559 * t676 + t670 * t748) * t577) * t733;
t538 = (((-t687 + t705) * t574 + t682 * t577) * t574 + t559 * t577 + (pkin(4) * t781 - t547) * t586) * t688;
t771 = 0.2e1 * t577;
t544 = (-pkin(4) * t586 + t670 * t771 + 0.2e1 * t781) * t586 * t628;
t631 = qJ(3,1) * t651 + t707;
t562 = t631 * t574;
t571 = t574 ^ 2;
t640 = qJ(3,1) * mrSges(3,2) + t762;
t619 = -t640 * t676 + t641 * t670;
t711 = qJ(3,1) * t774;
t775 = -m(3) * t687 + t711;
t622 = t692 + t775;
t722 = t670 * t676;
t644 = m(3) * qJ(3,1) + mrSges(3,3);
t745 = t577 * t644;
t759 = mrSges(3,2) * t670;
t768 = -0.2e1 * t631;
t712 = (-t622 * t658 + t722 * t768 + t708 + t775) * t544 + t619 * t535 - t538 * t759 + 0.4e1 * (t562 - t745 / 0.2e1) * t740 + (0.2e1 * (-t574 * t622 + t577 * t651) * t739 + (mrSges(3,2) * t771 - t574 * t641) * t574) * t676 - t571 * t640 * t670 - 0.2e1 * (t562 - t745) * t586;
t778 = t677 * t712;
t706 = -2 * t767 - Ifges(3,2) - Ifges(2,3);
t581 = t584 ^ 2;
t704 = (t617 * t542 + (-(t683 + t691) * m(3) + t709 + t706) * t533 + t651 * t536 + 0.2e1 * t572 * t747 + (t620 * t726 + t656 * t770 + t629) * t581) * t684;
t582 = t585 ^ 2;
t703 = (t618 * t543 + (-(t685 + t691) * m(3) + t710 + t706) * t534 + t651 * t537 + 0.2e1 * t573 * t746 + (t621 * t724 + t657 * t769 + t630) * t582) * t686;
t583 = t586 ^ 2;
t702 = (t619 * t544 + (-(t687 + t691) * m(3) + t711 + t706) * t535 + t651 * t538 + 0.2e1 * t574 * t745 + (t622 * t722 + t658 * t768 + t631) * t583) * t688;
t701 = (-m(3) * t536 + t533 * t651 - t542 * t761 - t569 * t642 + (t642 * t656 - t651 * t726 - t642) * t581) * t684;
t700 = (-m(3) * t537 + t534 * t651 - t543 * t760 - t570 * t643 + (t643 * t657 - t651 * t724 - t643) * t582) * t686;
t699 = (-m(3) * t538 + t535 * t651 - t544 * t759 - t571 * t644 + (t644 * t658 - t651 * t722 - t644) * t583) * t688;
t1 = [t603 * t699 + t601 * t700 + t599 * t701 + (t589 * t702 + t650 * t778) * t628 + (t588 * t703 + t649 * t779) * t627 + (t587 * t704 + t648 * t780) * t626; t604 * t699 + t602 * t700 + t600 * t701 + (t592 * t702 - t647 * t778) * t628 + (t591 * t703 - t646 * t779) * t627 + (t590 * t704 - t645 * t780) * t626; t696 * t699 + t697 * t700 + t698 * t701 + (-t712 * t671 + t702 * t736) * t628 + (-t713 * t669 + t703 * t737) * t627 + (-t714 * t667 + t704 * t738) * t626;];
taucX  = t1;
