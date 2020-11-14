% Calculate inertia matrix for parallel robot
% P3RRPRR8V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:02:51
% EndTime: 2020-08-06 21:02:52
% DurationCPUTime: 1.08s
% Computational Cost: add. (3183->232), mult. (4332->368), div. (249->9), fcn. (2445->35), ass. (0->145)
t624 = 2 * pkin(1);
t650 = 2 * mrSges(3,1);
t649 = 2 * mrSges(3,3);
t587 = cos(pkin(7));
t648 = 0.2e1 * t587;
t647 = m(3) * pkin(1);
t646 = m(3) * pkin(2);
t645 = mrSges(3,2) * pkin(1);
t582 = qJ(2,3) + pkin(7);
t595 = sin(qJ(2,3));
t609 = 0.2e1 * qJ(2,3);
t612 = pkin(3) ^ 2;
t613 = pkin(2) ^ 2;
t627 = 2 * pkin(2) * pkin(3);
t514 = sin(t609 + pkin(7)) * t627 + t613 * sin(t609) + t612 * sin(0.2e1 * t582) + (sin(t582) * pkin(3) + pkin(2) * t595) * t624;
t644 = t514 / 0.2e1;
t583 = qJ(2,2) + pkin(7);
t597 = sin(qJ(2,2));
t610 = 0.2e1 * qJ(2,2);
t515 = sin(t610 + pkin(7)) * t627 + t613 * sin(t610) + t612 * sin(0.2e1 * t583) + (sin(t583) * pkin(3) + pkin(2) * t597) * t624;
t643 = t515 / 0.2e1;
t584 = qJ(2,1) + pkin(7);
t599 = sin(qJ(2,1));
t611 = 0.2e1 * qJ(2,1);
t516 = sin(t611 + pkin(7)) * t627 + t613 * sin(t611) + t612 * sin(0.2e1 * t584) + (sin(t584) * pkin(3) + pkin(2) * t599) * t624;
t642 = t516 / 0.2e1;
t641 = pkin(3) * cos(t582);
t640 = pkin(3) * cos(t583);
t639 = pkin(3) * cos(t584);
t586 = sin(pkin(7));
t638 = pkin(3) * t586;
t594 = Ifges(3,2) - Ifges(3,1);
t637 = (-qJ(3,1) - pkin(5));
t636 = (-qJ(3,2) - pkin(5));
t635 = (-qJ(3,3) - pkin(5));
t634 = mrSges(3,2) * t586;
t633 = Ifges(3,4) * t586;
t632 = t586 * mrSges(3,1);
t631 = Ifges(2,5) + (-mrSges(2,1) - t646) * pkin(5);
t557 = t587 * pkin(3) + pkin(2);
t601 = cos(qJ(2,3));
t630 = 0.1e1 / (t557 * t601 - t595 * t638) * (t595 * t557 + t601 * t638);
t603 = cos(qJ(2,2));
t629 = 0.1e1 / (t557 * t603 - t597 * t638) * (t597 * t557 + t603 * t638);
t605 = cos(qJ(2,1));
t628 = 0.1e1 / (t557 * t605 - t599 * t638) * (t599 * t557 + t605 * t638);
t581 = t587 ^ 2;
t554 = t594 * t581;
t626 = (m(3) * t613) - 0.2e1 * pkin(2) * t634;
t625 = pkin(2) * t650;
t623 = -0.2e1 * pkin(1) * (mrSges(2,2) + t632);
t578 = -pkin(6) + t635;
t571 = 0.1e1 / t578;
t622 = t571 * t630;
t579 = -pkin(6) + t636;
t572 = 0.1e1 / t579;
t621 = t572 * t629;
t580 = -pkin(6) + t637;
t573 = 0.1e1 / t580;
t620 = t573 * t628;
t619 = -pkin(5) * mrSges(2,2) + Ifges(2,6);
t617 = -t634 + t646;
t618 = pkin(1) * t650 * t587 + (mrSges(2,1) + t617) * t624;
t616 = mrSges(3,2) * t587 + t632;
t615 = 0.4e1 * Ifges(3,4) * t581 + 0.2e1 * (-(pkin(2) * mrSges(3,2)) - t586 * t594) * t587 - 0.2e1 * pkin(2) * t632 + (2 * Ifges(2,4)) - 0.2e1 * Ifges(3,4);
t614 = Ifges(2,1) + Ifges(3,2) + Ifges(1,3) + (2 * (mrSges(2,3) + mrSges(3,3)) * pkin(5)) + ((m(2) + m(3)) * pkin(1) ^ 2) + (m(2) * pkin(5) ^ 2) - t554;
t606 = cos(qJ(1,1));
t604 = cos(qJ(1,2));
t602 = cos(qJ(1,3));
t600 = sin(qJ(1,1));
t598 = sin(qJ(1,2));
t596 = sin(qJ(1,3));
t590 = legFrame(1,3);
t589 = legFrame(2,3);
t588 = legFrame(3,3);
t576 = t605 * pkin(2);
t575 = t603 * pkin(2);
t574 = t601 * pkin(2);
t570 = cos(t590);
t569 = cos(t589);
t568 = cos(t588);
t567 = sin(t590);
t566 = sin(t589);
t565 = sin(t588);
t560 = t576 + pkin(1);
t559 = t575 + pkin(1);
t558 = t574 + pkin(1);
t552 = t637 * mrSges(3,1) + Ifges(3,5);
t551 = t636 * mrSges(3,1) + Ifges(3,5);
t550 = t635 * mrSges(3,1) + Ifges(3,5);
t549 = -t637 * mrSges(3,2) - Ifges(3,6);
t548 = -t636 * mrSges(3,2) - Ifges(3,6);
t547 = -t635 * mrSges(3,2) - Ifges(3,6);
t543 = 0.1e1 / (t576 + t639);
t542 = 0.1e1 / (t575 + t640);
t541 = 0.1e1 / (t574 + t641);
t540 = -mrSges(3,1) * t587 - t617;
t539 = t567 * t606 + t570 * t600;
t538 = t566 * t604 + t569 * t598;
t537 = t565 * t602 + t568 * t596;
t536 = -t567 * t600 + t570 * t606;
t535 = -t566 * t598 + t569 * t604;
t534 = -t565 * t596 + t568 * t602;
t533 = t587 * t625 + Ifges(2,3) + Ifges(3,3) + t626;
t532 = t560 * t606 - t600 * t580;
t531 = t559 * t604 - t598 * t579;
t530 = t558 * t602 - t596 * t578;
t529 = t600 * t560 + t606 * t580;
t528 = t598 * t559 + t604 * t579;
t527 = t596 * t558 + t602 * t578;
t520 = 0.2e1 * t554 + (t625 + 0.4e1 * t633) * t587 + Ifges(2,2) - Ifges(2,1) + t626 - t594;
t519 = t540 * t605 + t616 * t599 - t647;
t518 = t540 * t603 + t616 * t597 - t647;
t517 = t540 * t601 + t616 * t595 - t647;
t513 = t529 * t570 + t532 * t567 + t539 * t639;
t512 = t528 * t569 + t531 * t566 + t538 * t640;
t511 = t527 * t568 + t530 * t565 + t537 * t641;
t510 = -t529 * t567 + t532 * t570 + t536 * t639;
t509 = -t528 * t566 + t531 * t569 + t535 * t640;
t508 = -t527 * t565 + t530 * t568 + t534 * t641;
t507 = (t549 * t586 + t552 * t587 + ((-m(3) * qJ(3,1) - mrSges(3,3)) * pkin(2)) + t631) * t599 + t605 * (-t549 * t587 + t552 * t586 + t619);
t506 = (t548 * t586 + t551 * t587 + ((-m(3) * qJ(3,2) - mrSges(3,3)) * pkin(2)) + t631) * t597 + t603 * (-t548 * t587 + t551 * t586 + t619);
t505 = (t547 * t586 + t550 * t587 + ((-m(3) * qJ(3,3) - mrSges(3,3)) * pkin(2)) + t631) * t595 + t601 * (-t547 * t587 + t550 * t586 + t619);
t504 = (t543 * m(3) * t642 + t519 * t628) * t573;
t503 = (t542 * m(3) * t643 + t518 * t629) * t572;
t502 = (t541 * m(3) * t644 + t517 * t630) * t571;
t501 = (m(3) * t513 + t519 * t539) * t573;
t500 = (m(3) * t512 + t518 * t538) * t572;
t499 = (m(3) * t511 + t517 * t537) * t571;
t498 = (m(3) * t510 + t519 * t536) * t573;
t497 = (m(3) * t509 + t518 * t535) * t572;
t496 = (m(3) * t508 + t517 * t534) * t571;
t495 = (-t599 * t645 - t633) * t648 + t599 * t623 + (t637 ^ 2 * m(3)) + (qJ(3,1) * t649) + t614 + (t520 * t605 + t615 * t599 + t618) * t605;
t494 = (-t597 * t645 - t633) * t648 + t597 * t623 + (t636 ^ 2 * m(3)) + (qJ(3,2) * t649) + t614 + (t520 * t603 + t615 * t597 + t618) * t603;
t493 = (-t595 * t645 - t633) * t648 + t595 * t623 + (t635 ^ 2 * m(3)) + (qJ(3,3) * t649) + t614 + (t520 * t601 + t615 * t595 + t618) * t601;
t492 = (t495 * t539 + t513 * t519) * t573;
t491 = (t494 * t538 + t512 * t518) * t572;
t490 = (t493 * t537 + t511 * t517) * t571;
t489 = (t495 * t536 + t510 * t519) * t573;
t488 = (t494 * t535 + t509 * t518) * t572;
t487 = (t493 * t534 + t508 * t517) * t571;
t486 = -t495 * t620 + (t507 - t516 * t573 * t519 / 0.2e1) * t543;
t485 = -t494 * t621 + (t506 - t515 * t572 * t518 / 0.2e1) * t542;
t484 = -t493 * t622 + (t505 - t514 * t571 * t517 / 0.2e1) * t541;
t1 = [m(4) - (-t489 * t536 - t498 * t510) * t573 - (-t488 * t535 - t497 * t509) * t572 - (-t487 * t534 - t496 * t508) * t571, -(-t489 * t539 - t498 * t513) * t573 - (-t488 * t538 - t497 * t512) * t572 - (-t487 * t537 - t496 * t511) * t571, -(-t489 * t628 + (-t498 * t642 + t536 * t507) * t543) * t573 - (-t488 * t629 + (-t497 * t643 + t535 * t506) * t542) * t572 - (-t487 * t630 + (-t496 * t644 + t534 * t505) * t541) * t571; -(-t492 * t536 - t501 * t510) * t573 - (-t491 * t535 - t500 * t509) * t572 - (-t490 * t534 - t499 * t508) * t571, m(4) - (-t492 * t539 - t501 * t513) * t573 - (-t491 * t538 - t500 * t512) * t572 - (-t490 * t537 - t499 * t511) * t571, -(-t492 * t628 + (-t501 * t642 + t539 * t507) * t543) * t573 - (-t491 * t629 + (-t500 * t643 + t538 * t506) * t542) * t572 - (-t490 * t630 + (-t499 * t644 + t537 * t505) * t541) * t571; -(t486 * t536 - t504 * t510) * t573 - (t485 * t535 - t503 * t509) * t572 - (t484 * t534 - t502 * t508) * t571, -(t486 * t539 - t504 * t513) * t573 - (t485 * t538 - t503 * t512) * t572 - (t484 * t537 - t502 * t511) * t571, -t484 * t622 - t485 * t621 - t486 * t620 + m(4) + (t543 * t533 - (-t504 * t642 + t507 * t628) * t573) * t543 + (t542 * t533 - (-t503 * t643 + t506 * t629) * t572) * t542 + (t541 * t533 - (-t502 * t644 + t505 * t630) * t571) * t541;];
MX  = t1;
