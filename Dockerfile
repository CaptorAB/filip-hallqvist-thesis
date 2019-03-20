FROM trzeci/emscripten
COPY . /src
RUN node --version
RUN npm --version
RUN npm install
RUN npm run build
RUN npm test
RUN npm start