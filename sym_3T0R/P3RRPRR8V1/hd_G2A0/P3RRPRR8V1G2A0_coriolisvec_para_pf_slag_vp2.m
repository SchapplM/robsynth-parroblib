% Calculate vector of centrifugal and Coriolis load for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
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
% Datum: 2022-11-04 17:05
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:04:26
% EndTime: 2022-11-04 17:04:29
% DurationCPUTime: 2.71s
% Computational Cost: add. (13980->270), mult. (23631->461), div. (3051->12), fcn. (20190->35), ass. (0->218)
t630 = pkin(4) + qJ(3,3);
t639 = sin(qJ(1,3));
t645 = cos(qJ(1,3));
t629 = cos(pkin(5));
t735 = t629 * pkin(2);
t591 = pkin(1) + t735;
t644 = cos(qJ(2,3));
t628 = sin(pkin(5));
t638 = sin(qJ(2,3));
t718 = t628 * t638;
t668 = pkin(2) * t718 - t591 * t644;
t558 = -t645 * t630 - t639 * t668;
t739 = pkin(2) * t628;
t568 = t638 * t591 + t644 * t739;
t635 = legFrame(3,2);
t609 = sin(t635);
t612 = cos(t635);
t530 = t558 * t612 + t609 * t568;
t531 = -t558 * t609 + t612 * t568;
t615 = t644 * pkin(1);
t619 = qJ(2,3) + pkin(5);
t762 = t615 + pkin(2) * cos(t619);
t562 = t639 * t630 + t645 * t762;
t622 = 0.1e1 / t630;
t651 = xDP(3);
t652 = xDP(2);
t653 = xDP(1);
t521 = (t530 * t653 + t531 * t652 + t562 * t651) * t622;
t771 = 0.2e1 * t521;
t631 = pkin(4) + qJ(3,2);
t641 = sin(qJ(1,2));
t647 = cos(qJ(1,2));
t646 = cos(qJ(2,2));
t640 = sin(qJ(2,2));
t717 = t628 * t640;
t667 = pkin(2) * t717 - t591 * t646;
t559 = -t647 * t631 - t641 * t667;
t569 = t640 * t591 + t646 * t739;
t636 = legFrame(2,2);
t610 = sin(t636);
t613 = cos(t636);
t532 = t559 * t613 + t610 * t569;
t533 = -t559 * t610 + t613 * t569;
t616 = t646 * pkin(1);
t620 = qJ(2,2) + pkin(5);
t761 = t616 + pkin(2) * cos(t620);
t563 = t641 * t631 + t647 * t761;
t623 = 0.1e1 / t631;
t522 = (t532 * t653 + t533 * t652 + t563 * t651) * t623;
t770 = 0.2e1 * t522;
t632 = pkin(4) + qJ(3,1);
t643 = sin(qJ(1,1));
t649 = cos(qJ(1,1));
t648 = cos(qJ(2,1));
t642 = sin(qJ(2,1));
t716 = t628 * t642;
t666 = pkin(2) * t716 - t591 * t648;
t560 = -t649 * t632 - t643 * t666;
t570 = t642 * t591 + t648 * t739;
t637 = legFrame(1,2);
t611 = sin(t637);
t614 = cos(t637);
t534 = t560 * t614 + t611 * t570;
t535 = -t560 * t611 + t614 * t570;
t617 = t648 * pkin(1);
t621 = qJ(2,1) + pkin(5);
t760 = t617 + pkin(2) * cos(t621);
t564 = t643 * t632 + t649 * t760;
t624 = 0.1e1 / t632;
t523 = (t534 * t653 + t535 * t652 + t564 * t651) * t624;
t769 = 0.2e1 * t523;
t710 = 0.2e1 * pkin(1);
t747 = 0.2e1 * t629;
t768 = t747 / 0.2e1;
t766 = t629 * mrSges(3,2);
t574 = 0.1e1 / t762;
t548 = (t609 * t653 + t612 * t652) * t574;
t545 = t548 ^ 2;
t659 = pkin(2) ^ 2;
t660 = pkin(1) ^ 2;
t711 = -t659 - t660;
t585 = t735 * t710 - t711;
t765 = t545 * t585;
t575 = 0.1e1 / t761;
t549 = (t610 * t653 + t613 * t652) * t575;
t546 = t549 ^ 2;
t764 = t546 * t585;
t576 = 0.1e1 / t760;
t550 = (t611 * t653 + t614 * t652) * t576;
t547 = t550 ^ 2;
t763 = t547 * t585;
t731 = t628 * mrSges(3,1);
t759 = pkin(1) * t731 - Ifges(2,4) + Ifges(3,4);
t601 = mrSges(3,2) * qJ(3,1) - Ifges(3,6);
t604 = mrSges(3,1) * qJ(3,1) - Ifges(3,5);
t745 = m(3) * qJ(3,1);
t608 = mrSges(3,3) + t745;
t758 = -t608 * pkin(1) + t601 * t628 - t604 * t629 + Ifges(2,5);
t600 = mrSges(3,2) * qJ(3,2) - Ifges(3,6);
t603 = mrSges(3,1) * qJ(3,2) - Ifges(3,5);
t744 = m(3) * qJ(3,2);
t607 = mrSges(3,3) + t744;
t757 = -t607 * pkin(1) + t600 * t628 - t603 * t629 + Ifges(2,5);
t599 = mrSges(3,2) * qJ(3,3) - Ifges(3,6);
t602 = mrSges(3,1) * qJ(3,3) - Ifges(3,5);
t743 = m(3) * qJ(3,3);
t606 = mrSges(3,3) + t743;
t756 = -t606 * pkin(1) + t599 * t628 - t602 * t629 + Ifges(2,5);
t633 = Ifges(3,2) - Ifges(3,1);
t654 = m(3) * t660;
t734 = mrSges(3,2) * t628;
t572 = -0.2e1 * pkin(1) * t734 - Ifges(2,1) + Ifges(2,2) - t633 + t654;
t618 = t629 ^ 2;
t715 = t633 * t618;
t732 = Ifges(3,4) * t628;
t746 = pkin(1) * mrSges(3,1);
t551 = (0.4e1 * t732 + 0.2e1 * t746) * t629 + t572 + 0.2e1 * t715;
t755 = -0.2e1 * mrSges(3,3);
t699 = t609 * t739;
t702 = t612 * t739;
t721 = t609 * t591;
t724 = t591 * t612;
t536 = (-t639 * t721 + t702) * t644 + (t639 * t699 + t724) * t638;
t539 = (t639 * t724 + t699) * t644 + (-t639 * t702 + t721) * t638;
t565 = 0.1e1 / t668;
t518 = (t645 * t651 - (t536 * t652 + t539 * t653) * t565) * t622;
t754 = 0.2e1 * t518;
t698 = t610 * t739;
t701 = t613 * t739;
t720 = t610 * t591;
t723 = t591 * t613;
t537 = (-t641 * t720 + t701) * t646 + (t641 * t698 + t723) * t640;
t540 = (t641 * t723 + t698) * t646 + (-t641 * t701 + t720) * t640;
t566 = 0.1e1 / t667;
t519 = (t647 * t651 - (t537 * t652 + t540 * t653) * t566) * t623;
t753 = 0.2e1 * t519;
t697 = t611 * t739;
t700 = t614 * t739;
t719 = t611 * t591;
t722 = t591 * t614;
t538 = (-t643 * t719 + t700) * t648 + (t643 * t697 + t722) * t642;
t541 = (t643 * t722 + t697) * t648 + (-t643 * t700 + t719) * t642;
t567 = 0.1e1 / t666;
t520 = (t649 * t651 - (t538 * t652 + t541 * t653) * t567) * t624;
t752 = 0.2e1 * t520;
t751 = -0.2e1 * t521;
t750 = -0.2e1 * t522;
t749 = -0.2e1 * t523;
t584 = pkin(1) * mrSges(3,2) + t633 * t628;
t733 = Ifges(3,4) * t618;
t561 = t584 * t629 - 0.2e1 * t733 + t759;
t748 = 0.2e1 * t561;
t730 = t518 * t584;
t729 = t518 * t638;
t728 = t519 * t584;
t727 = t519 * t640;
t726 = t520 * t584;
t725 = t520 * t642;
t577 = pkin(1) * t638 + pkin(2) * sin(t619);
t656 = 0.2e1 * qJ(2,3);
t693 = pkin(2) * t710;
t500 = (-t765 + ((-t659 * cos(0.2e1 * t619) - t660 * cos(t656) + (-cos(t656 + pkin(5)) - t629) * t693 + t711) * t518 / 0.2e1 + t762 * t771 + (-t518 * t630 + 0.2e1 * t548 * t577) * t630) * t518) * t622;
t512 = (-0.1e1 / (t615 + (t629 * t644 - t718) * pkin(2)) * t765 + (t518 * t668 + t771) * t518) * t622;
t686 = -t602 * t628 + Ifges(2,6);
t542 = -t756 * t638 + (t599 * t629 - t686) * t644;
t655 = m(3) * pkin(1);
t672 = -mrSges(3,1) * t629 + t734;
t573 = -t655 + t672;
t583 = t731 + t766;
t552 = -t573 * t644 - t638 * t583;
t625 = t644 ^ 2;
t709 = 0.2e1 * t732;
t661 = t629 * t709 - Ifges(2,1) - Ifges(3,2) - Ifges(1,3) + t715;
t662 = 0.4e1 * (t709 + t746) * t629 + 0.2e1 * t572 + 0.4e1 * t715;
t689 = t545 * t574 * t577;
t692 = t518 * t561 * t625;
t696 = t518 * t733;
t714 = (t638 * t644 * t748 - t551 * t625 + (t755 - t743) * qJ(3,3) + t661) * t512 - t542 * t689 + t552 * t500 + t521 * t606 * t754 - t686 * t545 * t638 + (-t662 * t729 * t644 + t759 * t754 + t730 * t747 - 0.4e1 * t692 - 0.4e1 * t696 + (t599 * t638 * t768 + t756 * t644) * t548) * t548;
t578 = pkin(1) * t640 + pkin(2) * sin(t620);
t657 = 0.2e1 * qJ(2,2);
t501 = (-t764 + ((-t659 * cos(0.2e1 * t620) - t660 * cos(t657) + (-cos(pkin(5) + t657) - t629) * t693 + t711) * t519 / 0.2e1 + t761 * t770 + (-t519 * t631 + 0.2e1 * t549 * t578) * t631) * t519) * t623;
t513 = (-0.1e1 / (t616 + (t629 * t646 - t717) * pkin(2)) * t764 + (t519 * t667 + t770) * t519) * t623;
t685 = -t603 * t628 + Ifges(2,6);
t543 = -t757 * t640 + (t600 * t629 - t685) * t646;
t553 = -t573 * t646 - t640 * t583;
t626 = t646 ^ 2;
t688 = t546 * t575 * t578;
t691 = t519 * t561 * t626;
t695 = t519 * t733;
t713 = (t640 * t646 * t748 - t551 * t626 + (t755 - t744) * qJ(3,2) + t661) * t513 - t543 * t688 + t553 * t501 + t522 * t607 * t753 - t685 * t546 * t640 + (-t662 * t727 * t646 + t759 * t753 + t728 * t747 - 0.4e1 * t691 - 0.4e1 * t695 + (t600 * t640 * t768 + t757 * t646) * t549) * t549;
t579 = t642 * pkin(1) + pkin(2) * sin(t621);
t658 = 0.2e1 * qJ(2,1);
t502 = (-t763 + ((-t659 * cos(0.2e1 * t621) - t660 * cos(t658) + (-cos(pkin(5) + t658) - t629) * t693 + t711) * t520 / 0.2e1 + t760 * t769 + (-t520 * t632 + 0.2e1 * t550 * t579) * t632) * t520) * t624;
t514 = (-0.1e1 / (t617 + (t629 * t648 - t716) * pkin(2)) * t763 + (t520 * t666 + t769) * t520) * t624;
t684 = -t604 * t628 + Ifges(2,6);
t544 = -t758 * t642 + (t601 * t629 - t684) * t648;
t554 = -t573 * t648 - t642 * t583;
t627 = t648 ^ 2;
t687 = t547 * t576 * t579;
t690 = t520 * t561 * t627;
t694 = t520 * t733;
t712 = (t642 * t648 * t748 - t551 * t627 + (t755 - t745) * qJ(3,1) + t661) * t514 - t544 * t687 + t554 * t502 + t523 * t608 * t752 - t684 * t547 * t642 + (-t662 * t725 * t648 + t759 * t752 + t726 * t747 - 0.4e1 * t690 - 0.4e1 * t694 + (t601 * t642 * t768 + t758 * t648) * t550) * t550;
t708 = -0.2e1 * t731;
t706 = t638 * t751;
t705 = t640 * t750;
t704 = t642 * t749;
t683 = t714 * t565;
t682 = t713 * t566;
t681 = t712 * t567;
t680 = -m(3) * t500 + t552 * t512 + (-t518 * t606 / 0.2e1 + (-t573 * t638 + t583 * t644) * t548) * t754;
t679 = -m(3) * t501 + t553 * t513 + (-t519 * t607 / 0.2e1 + (-t573 * t640 + t583 * t646) * t549) * t753;
t678 = -m(3) * t502 + t554 * t514 + (-t520 * t608 / 0.2e1 + (-t573 * t642 + t583 * t648) * t550) * t752;
t571 = t672 * t710 - Ifges(2,3) - Ifges(3,3) - t654;
t589 = t655 - t734;
t671 = t574 * ((0.2e1 * t692 + (t521 * t708 + t551 * t729 + t751 * t766) * t644 + 0.2e1 * t696 + (mrSges(3,1) * t706 - t730) * t629 + t589 * t706 - t518 * t759) * t518 + t542 * t512 - t571 * t689);
t670 = t575 * ((0.2e1 * t691 + (t522 * t708 + t551 * t727 + t750 * t766) * t646 + 0.2e1 * t695 + (mrSges(3,1) * t705 - t728) * t629 + t589 * t705 - t519 * t759) * t519 + t543 * t513 - t571 * t688);
t669 = t576 * ((0.2e1 * t690 + (t523 * t708 + t551 * t725 + t749 * t766) * t648 + 0.2e1 * t694 + (mrSges(3,1) * t704 - t726) * t629 + t589 * t704 - t520 * t759) * t520 + t544 * t514 - t571 * t687);
t1 = [t611 * t669 + t610 * t670 + t609 * t671 + (t534 * t678 - t541 * t681) * t624 + (t532 * t679 - t540 * t682) * t623 + (t530 * t680 - t539 * t683) * t622; t614 * t669 + t613 * t670 + t612 * t671 + (t535 * t678 - t538 * t681) * t624 + (t533 * t679 - t537 * t682) * t623 + (t531 * t680 - t536 * t683) * t622; (t564 * t678 + t649 * t712) * t624 + (t563 * t679 + t647 * t713) * t623 + (t562 * t680 + t645 * t714) * t622;];
taucX  = t1;
