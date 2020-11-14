% Calculate Gravitation load for parallel robot
% P3RPRRR12V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR12V1G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:27
% EndTime: 2020-08-06 18:21:27
% DurationCPUTime: 0.40s
% Computational Cost: add. (339->83), mult. (441->133), div. (21->7), fcn. (321->18), ass. (0->62)
t587 = m(2) + m(3);
t586 = pkin(1) + pkin(5);
t601 = -mrSges(1,2) + mrSges(2,3);
t569 = legFrame(3,3);
t559 = sin(t569);
t562 = cos(t569);
t542 = -t559 * g(1) + t562 * g(2);
t545 = t562 * g(1) + t559 * g(2);
t548 = m(2) * pkin(1) + t586 * m(3) + mrSges(1,1) - mrSges(2,2) + mrSges(3,3);
t572 = sin(qJ(3,3));
t556 = t572 * pkin(3) + qJ(2,3);
t553 = 0.1e1 / t556;
t573 = sin(qJ(1,3));
t579 = cos(qJ(1,3));
t578 = cos(qJ(3,3));
t591 = mrSges(3,1) * t572 + mrSges(3,2) * t578 + t587 * qJ(2,3) + t601;
t600 = ((t548 * t573 - t591 * t579) * t545 + (-t548 * t579 - t591 * t573) * t542) * t553;
t570 = legFrame(2,3);
t560 = sin(t570);
t563 = cos(t570);
t543 = -t560 * g(1) + t563 * g(2);
t546 = t563 * g(1) + t560 * g(2);
t574 = sin(qJ(3,2));
t557 = t574 * pkin(3) + qJ(2,2);
t554 = 0.1e1 / t557;
t575 = sin(qJ(1,2));
t581 = cos(qJ(1,2));
t580 = cos(qJ(3,2));
t590 = mrSges(3,1) * t574 + mrSges(3,2) * t580 + t587 * qJ(2,2) + t601;
t599 = ((t548 * t575 - t590 * t581) * t546 + (-t548 * t581 - t590 * t575) * t543) * t554;
t571 = legFrame(1,3);
t561 = sin(t571);
t564 = cos(t571);
t544 = -t561 * g(1) + t564 * g(2);
t547 = t564 * g(1) + t561 * g(2);
t576 = sin(qJ(3,1));
t558 = t576 * pkin(3) + qJ(2,1);
t555 = 0.1e1 / t558;
t577 = sin(qJ(1,1));
t583 = cos(qJ(1,1));
t582 = cos(qJ(3,1));
t589 = mrSges(3,1) * t576 + mrSges(3,2) * t582 + t587 * qJ(2,1) + t601;
t598 = ((t548 * t577 - t589 * t583) * t547 + (-t548 * t583 - t589 * t577) * t544) * t555;
t594 = t579 * t542 - t545 * t573;
t597 = t594 * t553;
t593 = t581 * t543 - t546 * t575;
t596 = t593 * t554;
t592 = t583 * t544 - t547 * t577;
t595 = t592 * t555;
t585 = mrSges(3,1) * g(3);
t584 = mrSges(3,2) * g(3);
t568 = 0.1e1 / t576;
t567 = 0.1e1 / t574;
t566 = 0.1e1 / t572;
t565 = pkin(6) + t586;
t541 = -t558 * t583 + t565 * t577;
t540 = -t557 * t581 + t565 * t575;
t539 = -t556 * t579 + t565 * t573;
t538 = t577 * t558 + t565 * t583;
t537 = t575 * t557 + t565 * t581;
t536 = t573 * t556 + t565 * t579;
t1 = [(-t561 * t577 + t564 * t583) * t598 + (-t560 * t575 + t563 * t581) * t599 + (-t559 * t573 + t562 * t579) * t600 - g(1) * m(4) + ((t538 * t564 - t561 * t541) * t595 + (t537 * t563 - t560 * t540) * t596 + (t536 * t562 - t559 * t539) * t597) * t587; (t561 * t583 + t564 * t577) * t598 + (t560 * t581 + t563 * t575) * t599 + (t559 * t579 + t562 * t573) * t600 - g(2) * m(4) + ((t561 * t538 + t541 * t564) * t595 + (t560 * t537 + t540 * t563) * t596 + (t559 * t536 + t539 * t562) * t597) * t587; -g(3) * m(4) + (-t568 * ((t592 * mrSges(3,1) + t584) * t582 + t576 * (-t592 * mrSges(3,2) + t585)) - t567 * ((t593 * mrSges(3,1) + t584) * t580 + t574 * (-t593 * mrSges(3,2) + t585)) - t566 * ((t594 * mrSges(3,1) + t584) * t578 + t572 * (-t594 * mrSges(3,2) + t585))) / pkin(3) + (t566 * t594 * t578 + t567 * t593 * t580 + t568 * t592 * t582) * t587;];
taugX  = t1;
