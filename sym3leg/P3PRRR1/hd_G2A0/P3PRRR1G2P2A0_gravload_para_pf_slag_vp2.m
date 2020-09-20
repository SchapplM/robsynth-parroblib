% Calculate Gravitation load for parallel robot
% P3PRRR1G2P2A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR1G2P2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:23
% EndTime: 2020-03-09 21:18:23
% DurationCPUTime: 0.22s
% Computational Cost: add. (474->68), mult. (399->106), div. (45->5), fcn. (291->18), ass. (0->60)
t601 = pkin(7) + qJ(2,3);
t591 = qJ(3,3) + t601;
t578 = sin(t591);
t605 = legFrame(3,2);
t594 = sin(t605);
t597 = cos(t605);
t575 = t597 * g(1) - t594 * g(2);
t581 = cos(t591);
t608 = mrSges(3,2) * g(3);
t610 = mrSges(3,1) * g(3);
t551 = -t581 * (mrSges(3,1) * t575 - t608) + (mrSges(3,2) * t575 + t610) * t578;
t600 = m(3) * pkin(2) + mrSges(2,1);
t584 = t600 * g(3);
t585 = sin(t601);
t588 = cos(t601);
t609 = mrSges(2,2) * g(3);
t624 = 0.1e1 / (-t588 * t578 + t581 * t585);
t618 = t624 * ((-t575 * t600 + t609) * t588 + t585 * (t575 * mrSges(2,2) + t584) + t551);
t627 = t578 * t618;
t602 = pkin(7) + qJ(2,2);
t592 = qJ(3,2) + t602;
t579 = sin(t592);
t606 = legFrame(2,2);
t595 = sin(t606);
t598 = cos(t606);
t576 = t598 * g(1) - t595 * g(2);
t582 = cos(t592);
t552 = -t582 * (mrSges(3,1) * t576 - t608) + (mrSges(3,2) * t576 + t610) * t579;
t586 = sin(t602);
t589 = cos(t602);
t623 = 0.1e1 / (-t589 * t579 + t582 * t586);
t617 = t623 * ((-t576 * t600 + t609) * t589 + t586 * (t576 * mrSges(2,2) + t584) + t552);
t626 = t579 * t617;
t603 = pkin(7) + qJ(2,1);
t593 = qJ(3,1) + t603;
t580 = sin(t593);
t607 = legFrame(1,2);
t596 = sin(t607);
t599 = cos(t607);
t577 = t599 * g(1) - t596 * g(2);
t583 = cos(t593);
t553 = -t583 * (mrSges(3,1) * t577 - t608) + (mrSges(3,2) * t577 + t610) * t580;
t587 = sin(t603);
t590 = cos(t603);
t622 = 0.1e1 / (-t590 * t580 + t583 * t587);
t616 = t622 * ((-t577 * t600 + t609) * t590 + t587 * (t577 * mrSges(2,2) + t584) + t553);
t625 = t580 * t616;
t621 = t551 * t624;
t620 = t552 * t623;
t619 = t553 * t622;
t615 = (pkin(2) * t585 + pkin(3) * t578) * t621;
t614 = (pkin(2) * t586 + pkin(3) * t579) * t620;
t613 = (pkin(2) * t587 + pkin(3) * t580) * t619;
t612 = 0.1e1 / pkin(2);
t611 = 0.1e1 / pkin(3);
t604 = m(1) + m(2) + m(3);
t574 = t596 * g(1) + t599 * g(2);
t573 = t595 * g(1) + t598 * g(2);
t572 = t594 * g(1) + t597 * g(2);
t1 = [-g(1) * m(4) + (-t594 * t572 - t595 * t573 - t596 * t574) * t604 + (-t597 * t627 - t598 * t626 - t599 * t625 + (t597 * t615 + t598 * t614 + t599 * t613) * t611) * t612; -g(2) * m(4) + (-t572 * t597 - t573 * t598 - t574 * t599) * t604 + (t594 * t627 + t595 * t626 + t596 * t625 + (-t594 * t615 - t595 * t614 - t596 * t613) * t611) * t612; -g(3) * m(4) + (-t581 * t618 - t582 * t617 - t583 * t616 + ((pkin(2) * t590 + pkin(3) * t583) * t619 + (pkin(2) * t589 + pkin(3) * t582) * t620 + (pkin(2) * t588 + pkin(3) * t581) * t621) * t611) * t612;];
taugX  = t1;
